clear;

Nr = 3;
rr = 2;

N = 1000;

A = zeros(N,N);
K = zeros(N,N);

G = Ggen_alter(rr,N,Nr);

for i = 1:N
    for j = 1:N
      A(i,j) = (abs(G(i,j)))^2; 
    end
end

K = A - eye(N).*sum(A,2); 

evals = eig(K);

figure(1)
plot(real(evals),imag(evals),'ro','MarkerSize',3)
xlim([-1.200 -0.8])
xlabel('Re(\lambda)', 'FontName', 'Times New Roman','FontSize',20);
ylabel('Im(\lambda)', 'FontName', 'Times New Roman','FontSize',20);
set(gca,'FontSize',16)
set(gca,'LineWidth',1.5);
