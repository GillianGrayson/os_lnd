%str1 = getenv('SLURM_ARRAY_TASK_ID'); %raskomentirovat dlja klastera (takzhe kak i poslednyuu stroku)
str1 = '1'; %zakomentirovat dlja klastera
num = str2num(str1);

%%%% parameters %%%%

Nr = 1; % number of realizations of a random state
N = 50;
M = N^2; % size of rho-matrix
p = 0.1; % parameter for partial decoherence 0<=p<=1

rho_temp = St_gen(num,M,Nr); % generation of random state 

rho = Decoh(rho_temp,M,p); % decoherence of the random state

%******************************************************




%**********************************************************
F = Reshuff(rho,N,M); % reshuffling procedure

%%% ****************8calculation of F^*(id)
y1 = zeros(M,1);
y2 = eye(N);
%y1 = reshape(y2,[],1);

S2 = zeros(N,N);
for s1 = 1:N
for s2 = 1:N
kk=s2+(N-1)*s1;
y1(kk)=y2(s1,s2);
end
end


FT=ctranspose(F);
S=FT*y1;
%*************************8
%S1=vec2mat(S,N);

S2 = zeros(N,N);
for s1 = 1:N
for s2 = 1:N
kk=s2+(N-1)*s1;
S2(s1,s2)=S(kk);
end
end

trace (S2)

%size(S1)
%TTT= kron(S1,eye(N));
%size(TTT)


%%%%% generation of Lindblad operator L %%%%%%
L = F - (1/2)*(kron(S2,eye(N)) + kron(eye(N),transpose(S2))); 

evals = eig(L); % calculation of Lindblad eigenvalues 

Evals_re = real(evals);
Evals_im = imag(evals);

figure
scatter(Evals_re(2:end), Evals_im(2:end))

hull = fopen('sobstvenie_znacheniya.dat','w');
%fprintf(hull,'%8.4f %16.12f\n',DT.Points(C,1),DT.Points(C,2));
for nn=1:M
sx=Evals_re(nn);
sy=Evals_im(nn);
fprintf(hull,'%12.8f %16.12f\n',sx,sy);
end


fclose(hull);


%quit() %raskomentirovat dlja klastera