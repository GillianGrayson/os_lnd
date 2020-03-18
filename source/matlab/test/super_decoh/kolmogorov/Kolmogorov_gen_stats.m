str1 = getenv('SLURM_ARRAY_TASK_ID');
%str1 = '1';
number = str2num(str1);

Nr = 1000;
N = 1000;

A = zeros(N,N);
K = zeros(N,N);

G = Ggen_alter(number,N,Nr);

for i = 1:N
    for j = 1:N
      A(i,j) = (abs(G(i,j)))^2; 
    end
end

K = A - eye(N).*sum(A,2); 

evals = eig(K);

Evals_re = real(evals);
Evals_im = imag(evals);

fname1 = strcat('evals_K_re_',str1,'.txt');
fname2 = strcat('evals_K_im_',str1,'.txt');

%write evals of L-matrix to txt-files
save(fname1, 'Evals_re', '-ASCII', '-double')
save(fname2, 'Evals_im', '-ASCII', '-double')

quit()