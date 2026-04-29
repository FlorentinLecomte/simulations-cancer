% Linear stability analysis of the model of the article
% Plot of Lambda_k as a function of n
% If negative, a_k->0, otherwise a_k grows exponentially

% Parameters
dn=0.001;
dm=0.001;
dn_mot=0.01; %increased motility
dm_mot=0.01;
gamma=0.005;
eta=10;
alpha=0.1;
beta=1;
epsilon=0.01;
mu=10;
nu=10;

% Other parameters
n=1;
m=alpha/beta;
k=1:100;
lam_k=zeros(size(k));
lam_k_mot=zeros(size(k));

for i=1:length(k)
    % Jacobian matrices
    Jr= [-mu, -mu, 0;
        -nu, -eta*m-nu, 0;
        alpha, 0,-beta];
    Jt= [dn, -gamma*n, 0;
        0, 0, 0;
        0, 0, dm];

    % Increased motility
    Jt_mot= [dn_mot, -gamma*n, 0;
        0, 0, 0;
        0, 0, dm_mot];

    % Lambda_k
    eigenvalues= eig(Jr-k(i)^2*Jt);
    eigenvalues_mot= eig(Jr-k(i)^2*Jt_mot);
    lam_k(i)=max(real(eigenvalues));
    lam_k_mot(i)=max(real(eigenvalues_mot));
end

% Plot
figure;
plot(k, lam_k, 'r');
hold on;
plot(k, lam_k_mot, 'b');
ylim([-10,0]);
title('Stability as a function of k with n=1 in 1D');
xlabel('k values');
ylabel('\lambda_k');
legend('Initial parameters', 'Increased motility', 'Location', 'northeast');
hold off;

saveas(gcf, 'figures/lin_stab_analysis2.png');