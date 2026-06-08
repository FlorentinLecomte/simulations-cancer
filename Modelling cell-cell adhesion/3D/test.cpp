#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <omp.h>
#include <algorithm>

using namespace std;

// --- Paramètres du papier (Table 1 & 3) ---
const double R_MAX = 5.0;
const double E_MODULUS = 1e-3;
const double POISSON = 0.5;
const double ALPHA_ADH = 3.722e-4;
const double GAMMA_FRIC = 0.01;
const double A_T = 4e-3;
const double T_CC = 1000.0;
const double DT = 1.0;
const double T_MAX = 4000.0;
const double T_SAVE = 50.0;

// Paramètres des fibres (Table 3)
const double ALPHA_FIBRE = 0.01;     // Coefficient d'adhésion cellule-fibre
const double BETA_FIBRE = 1e-3;      // Coefficient résistif
const double V_MAX = 10.0;           // Vitesse max induite par la fibre (µm/min)
const double FIBRE_RADIUS = 2.0;     // Rayon de la fibre (µm)
const double P_CONTACT = 1e-3;       // Probabilité de dégradation (ajustée pour la stabilité)
const double P_DEGRADATION = 0.05; // Probabilité par seconde de dégrader une fibre touchée
const double P_DEPOSITION = 0.02;  // Probabilité par seconde de déposer une fibre
const int MAX_FIBRES = 50000;     // Taille de notre Object Pool

struct Cell {
    double x, y, z;
    double vx, vy, vz;       // On doit maintenant conserver la vitesse pour les fibres
    double fx, fy, fz;
    double radius;
    int num_neighbours;
    bool is_dead;
};

struct Fibre {
    double x1, y1, z1;       // Point de départ
    double x2, y2, z2;       // Point d'arrivée
    double dx, dy, dz;       // Vecteur directionnel normalisé (l_f)
    double length;
    bool is_degraded;
};

double calculate_effective_E() {
    double term = (1.0 - (POISSON * POISSON)) / E_MODULUS;
    return 1.0 / (term + term);
}

double calculate_effective_R(double r1, double r2) {
    return (r1 * r2) / (r1 + r2);
}

// Fonction pour calculer la distance entre une cellule et un segment de fibre
double distance_point_segment(double px, double py, double pz, const Fibre& f, double& proj_x, double& proj_y, double& proj_z) {
    double v_x = px - f.x1;
    double v_y = py - f.y1;
    double v_z = pz - f.z1;

    double u_x = f.x2 - f.x1;
    double u_y = f.y2 - f.y1;
    double u_z = f.z2 - f.z1;

    double dot_vu = (v_x * u_x + v_y * u_y + v_z * u_z);
    double dot_uu = (u_x * u_x + u_y * u_y + u_z * u_z);

    double t = dot_vu / dot_uu;
    t = std::max(0.0, std::min(1.0, t)); // Restreindre au segment

    proj_x = f.x1 + t * u_x;
    proj_y = f.y1 + t * u_y;
    proj_z = f.z1 + t * u_z;

    double dist_x = px - proj_x;
    double dist_y = py - proj_y;
    double dist_z = pz - proj_z;

    return std::sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
}

int main() {
    mt19937 gen(42);
    uniform_real_distribution<double> dist_pos(-50.0, 50.0); // Domaine plus grand
    uniform_real_distribution<double> prob_dist(0.0, 1.0);
    normal_distribution<double> norm_dist(0.0, 1.0);

    vector<Cell> cells;
    cells.reserve(50000);

    // Initialisation amas central
    for(int i = 0; i < 10; i++) {
        cells.push_back({
            dist_pos(gen)*0.1, dist_pos(gen)*0.1, dist_pos(gen)*0.1,
            0,0,0, 0,0,0, R_MAX * 0.9, 0, false
        });
    }

    // Initialisation des fibres (alignées sur l'axe Y)
    int num_fibres_initiales = 2000;
    vector<Fibre> fibres;
    fibres.reserve(MAX_FIBRES); // Optimisation de la mémoire

    // 1. Créer les 2000 fibres initiales actives
    for(int i = 0; i < num_fibres_initiales; i++) {
        double fx = dist_pos(gen) * 2.0;
        double fy = dist_pos(gen) * 2.0;
        double fz = dist_pos(gen) * 2.0;
        double length = 75.0 + norm_dist(gen) * 5.0;

        fibres.push_back({
            fx, fy - length/2.0, fz,
            fx, fy + length/2.0, fz,
            0, 1.0, 0,
            length, false // false = active
        });
    }

    // 2. Remplir le reste du pool avec des fibres inactives
    for(int i = num_fibres_initiales; i < MAX_FIBRES; i++) {
        fibres.push_back({
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            0, true // true = dégradée (inactive et prête à être recyclée)
        });
    }

    ofstream f_data("fibres.bin", ios::binary);
    if (f_data.is_open()) {
        for(int i = 0; i < MAX_FIBRES; i++) {
            f_data.write(reinterpret_cast<const char*>(&fibres[i].x1), sizeof(double));
            f_data.write(reinterpret_cast<const char*>(&fibres[i].y1), sizeof(double));
            f_data.write(reinterpret_cast<const char*>(&fibres[i].z1), sizeof(double));
            f_data.write(reinterpret_cast<const char*>(&fibres[i].x2), sizeof(double));
            f_data.write(reinterpret_cast<const char*>(&fibres[i].y2), sizeof(double));
            f_data.write(reinterpret_cast<const char*>(&fibres[i].z2), sizeof(double));
        }
        f_data.close();
    }

    double E_star = calculate_effective_E();

    ofstream data("data.bin", ios::binary);
    if (!data.is_open()) return 1;

    int save_interval = round(T_SAVE / DT);
    int current_step = 0;
    double t1 = omp_get_wtime();

    // --- BOUCLE DE SIMULATION ---
    for(double t = 0; t <= T_MAX; t += DT) {
        int n_cells = cells.size();

        #pragma omp parallel for
        for(int i = 0; i < n_cells; i++) {
            cells[i].fx = 0; cells[i].fy = 0; cells[i].fz = 0;
            cells[i].num_neighbours = 0;
        }

        // --- 1. Forces Cellule-Cellule ---
        #pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < n_cells; i++) {
            if(cells[i].is_dead) continue;

            for(int j = i + 1; j < n_cells; j++) {
                if(cells[j].is_dead) continue;

                double dx = cells[j].x - cells[i].x;
                double dy = cells[j].y - cells[i].y;
                double dz = cells[j].z - cells[i].z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                if(dist < (cells[i].radius + cells[j].radius)) {
                    double h_ij = cells[i].radius + cells[j].radius - dist;
                    double R_star = calculate_effective_R(cells[i].radius, cells[j].radius);

                    double F_rep = (4.0/3.0) * E_star * std::sqrt(R_star) * std::pow(h_ij, 1.5);
                    double F_adh = ALPHA_ADH * (cells[i].radius - (h_ij / 4.0)) * h_ij;
                    double F_total_mag = F_rep - F_adh;

                    double nx = dx / dist, ny = dy / dist, nz = dz / dist;

                    #pragma omp atomic
                    cells[i].fx -= F_total_mag * nx;
                    #pragma omp atomic
                    cells[i].fy -= F_total_mag * ny;
                    #pragma omp atomic
                    cells[i].fz -= F_total_mag * nz;
                    #pragma omp atomic
                    cells[i].num_neighbours++;

                    #pragma omp atomic
                    cells[j].fx += F_total_mag * nx;
                    #pragma omp atomic
                    cells[j].fy += F_total_mag * ny;
                    #pragma omp atomic
                    cells[j].fz += F_total_mag * nz;
                    #pragma omp atomic
                    cells[j].num_neighbours++;
                }
            }
        }

        // --- 2. Forces Cellule-Fibre ---
        // Ajoute ceci au début de ta boucle de simulation ou gère un générateur par thread
        #pragma omp parallel
        {
            // Chaque thread a son propre générateur pour éviter les goulets d'étranglement
            mt19937 local_gen(42 + omp_get_thread_num());
            uniform_real_distribution<double> local_prob(0.0, 1.0);

            #pragma omp for schedule(dynamic)
            for(int i = 0; i < n_cells; i++) {
                if(cells[i].is_dead) continue;

                for(int f = 0; f < MAX_FIBRES; f++) { // On boucle sur le pool entier
                    if(fibres[f].is_degraded) continue;

                    double proj_x, proj_y, proj_z;
                    double dist = distance_point_segment(cells[i].x, cells[i].y, cells[i].z, fibres[f], proj_x, proj_y, proj_z);

                    if(dist < (cells[i].radius + FIBRE_RADIUS)) {
                        // Vitesse parallèle à la fibre (v_i . l_f)
                        double v_dot_l = cells[i].vx * fibres[f].dx + cells[i].vy * fibres[f].dy + cells[i].vz * fibres[f].dz;

                        // Calcul de la composante d'adhésion parallèle (Eq 15)
                        double adhesion_factor = ALPHA_FIBRE * (1.0 - (v_dot_l / V_MAX));
                        double f_para_x = adhesion_factor * fibres[f].dx;
                        double f_para_y = adhesion_factor * fibres[f].dy;
                        double f_para_z = adhesion_factor * fibres[f].dz;

                        // Calcul de la vélocité orthogonale pour la friction (Eq 16)
                        double v_ortho_x = cells[i].vx - v_dot_l * fibres[f].dx;
                        double v_ortho_y = cells[i].vy - v_dot_l * fibres[f].dy;
                        double v_ortho_z = cells[i].vz - v_dot_l * fibres[f].dz;

                        double f_perp_x = BETA_FIBRE * v_ortho_x;
                        double f_perp_y = BETA_FIBRE * v_ortho_y;
                        double f_perp_z = BETA_FIBRE * v_ortho_z;

                        cells[i].fx += (f_para_x - f_perp_x);
                        cells[i].fy += (f_para_y - f_perp_y);
                        cells[i].fz += (f_para_z - f_perp_z);

                        // NOUVEAU : Logique de dégradation enzymatique (MT1-MMP)
                        if (local_prob(local_gen) < P_DEGRADATION * DT) {
                            fibres[f].is_degraded = true;
                            // Note: En multithreading strict, écrire un booléen à true depuis
                            // plusieurs threads est bénin, pas besoin de #pragma omp atomic ici.
                        }
                    }
                }
            }
        }

        // --- 3. Intégration du mouvement et Croissance ---
        for(int i = 0; i < n_cells; i++) {
            if(cells[i].is_dead) continue;

            double rand_fx = A_T * norm_dist(gen);
            double rand_fy = A_T * norm_dist(gen);
            double rand_fz = A_T * norm_dist(gen);

            // Calcul et stockage de la nouvelle vitesse
            cells[i].vx = (cells[i].fx + rand_fx) / GAMMA_FRIC;
            cells[i].vy = (cells[i].fy + rand_fy) / GAMMA_FRIC;
            cells[i].vz = (cells[i].fz + rand_fz) / GAMMA_FRIC;

            // Mise à jour de la position
            cells[i].x += cells[i].vx * DT;
            cells[i].y += cells[i].vy * DT;
            cells[i].z += cells[i].vz * DT;

            // NOUVEAU : Logique de dépôt de matrice
            if(prob_dist(gen) < P_DEPOSITION * DT) { // (Ici on est hors du #pragma omp, on utilise gen global)

                // Normaliser la vitesse pour orienter la nouvelle fibre
                double v_mag = std::sqrt(cells[i].vx*cells[i].vx + cells[i].vy*cells[i].vy + cells[i].vz*cells[i].vz);
                if (v_mag > 0.001) {
                    double dir_x = cells[i].vx / v_mag;
                    double dir_y = cells[i].vy / v_mag;
                    double dir_z = cells[i].vz / v_mag;

                    // Chercher une fibre "morte" dans le pool pour la recycler (Object Pooling)
                    for(int f = 0; f < MAX_FIBRES; f++) {
                        if(fibres[f].is_degraded) {
                            double new_length = 30.0; // Les fibres déposées sont souvent plus courtes
                            fibres[f].x1 = cells[i].x - dir_x * (new_length/2.0);
                            fibres[f].y1 = cells[i].y - dir_y * (new_length/2.0);
                            fibres[f].z1 = cells[i].z - dir_z * (new_length/2.0);

                            fibres[f].x2 = cells[i].x + dir_x * (new_length/2.0);
                            fibres[f].y2 = cells[i].y + dir_y * (new_length/2.0);
                            fibres[f].z2 = cells[i].z + dir_z * (new_length/2.0);

                            fibres[f].dx = dir_x;
                            fibres[f].dy = dir_y;
                            fibres[f].dz = dir_z;
                            fibres[f].length = new_length;
                            fibres[f].is_degraded = false; // La fibre "revient à la vie"

                            break; // On a trouvé une place, on sort de la boucle
                        }
                    }
                }
            }

            if(cells[i].radius < R_MAX) {
                cells[i].radius += 0.1 * DT;
                if(cells[i].radius > R_MAX) cells[i].radius = R_MAX;
            }
        }

        // --- 4. Mitose ---
        vector<Cell> new_cells;
        for(int i = 0; i < n_cells; i++) {
            if(cells[i].is_dead) continue;
            if(cells[i].radius >= 0.99 * R_MAX && cells[i].num_neighbours < 16) {
                if(prob_dist(gen) < (1.0 / T_CC) * DT) {
                    cells[i].radius = R_MAX * 0.8;
                    Cell daughter = cells[i];
                    daughter.x += 0.1 * norm_dist(gen);
                    daughter.y += 0.1 * norm_dist(gen);
                    daughter.z += 0.1 * norm_dist(gen);
                    daughter.num_neighbours = 0;
                    new_cells.push_back(daughter);
                }
            }
        }
        cells.insert(cells.end(), new_cells.begin(), new_cells.end());

        // --- 5. Sauvegarde ---
        if (current_step % save_interval == 0) {
            cout << "Simulation à t = " << t << " min | Population: " << cells.size() << " cellules" << endl;
            for(int i = 0; i < cells.size(); i++) {
                if(cells[i].is_dead) continue;
                double val_to_plot = cells[i].radius;
                data.write(reinterpret_cast<const char*>(&t), sizeof(double));
                data.write(reinterpret_cast<const char*>(&cells[i].x), sizeof(double));
                data.write(reinterpret_cast<const char*>(&cells[i].y), sizeof(double));
                data.write(reinterpret_cast<const char*>(&cells[i].z), sizeof(double));
                data.write(reinterpret_cast<const char*>(&val_to_plot), sizeof(double));
            }
        }
        current_step++;
    }

    data.close();
    cout << "Simulation 3D terminée en " << (omp_get_wtime() - t1) << " secondes." << endl;
    return 0;
}