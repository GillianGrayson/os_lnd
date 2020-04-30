%str1 = getenv('SLURM_ARRAY_TASK_ID'); %raskomentirovat dlja klastera (takzhe kak i poslednyuu stroku)
str1 = '1'; %zakomentirovat dlja klastera
num = str2num(str1);

%%%% parameters %%%%

Nr = 1; % number of realizations of a random state
N = 10;
M = N^2; % size of rho-matrix
p = 0.8; % parameter for partial decoherence 0<=p<=1

hull = fopen('sobstvenie_znacheniya.dat','w');

imag1 = sqrt(-1.);

for jjj=1:100

%rho_temp = St_gen(num,M,Nr); % generation of random state 

X=zeros(M);

for s = 1:M
for sp = 1:M
temp_y = randn;
temp_z = randn;
X(s,sp) = (temp_y + 1.*imag1*temp_z)/sqrt(2);
end
end
rho_temp=X*X';
rho_temp = sqrt(M)*rho_temp/trace(rho_temp);


rho = Decoh(rho_temp,M,p); % decoherence of the random state

%******************************************************




%**********************************************************
F = Reshuff(rho,N,M); % reshuffling procedure

AA=zeros(N);

for s1 = 1:N
for s2 = 1:N
for s3 = 1:N
w1=s3+N*(s1-1);
w2=s3+N*(s2-1);
AA(s1,s2)=AA(s1,s2)+rho(w1,w2);
end
end
end



%%%%% generation of Lindblad operator L %%%%%%
L = F - (1/2)*(kron(AA,eye(N)) + kron(eye(N),transpose(AA))); 

evals = eig(L); % calculation of Lindblad eigenvalues 

Evals_re = real(evals);
Evals_im = imag(evals);




%fprintf(hull,'%8.4f %16.12f\n',DT.Points(C,1),DT.Points(C,2));
for nn=2:M
sx=N*(Evals_re(nn)+1);
sy=N*Evals_im(nn);
fprintf(hull,'%12.8f %16.12f\n',sx,sy);
end

jjj

end
fclose(hull);


%quit() %raskomentirovat dlja klastera
