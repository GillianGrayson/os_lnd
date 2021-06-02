%str1 = getenv('SLURM_ARRAY_TASK_ID'); %raskomentirovat dlja klastera (takzhe kak i poslednyuu stroku)
str1 = '1'; %zakomentirovat dlja klastera
num = str2num(str1);

%%%% parameters %%%%

Nr = 1; % number of realizations of a random state
N = 10;
M = N^2; % size of rho-matrix
hull = fopen('sobstvenie_znacheniya_Veleva80.dat','w');

imag1 = sqrt(-1.);


%*****************************************************************************************


%%% calculate MR-matrix
seed = num;
rand('state', seed);
MR = randn(N,'double');

%%% calculate MI-matrix
MI = randn(N,'double');

MS = zeros(N,N);


MS=MR+imag1*MI;

MZ = zeros(N,N);

for j = 1:N
    for k = 1:N
        MZ(k,j) = MS(k,j)*conj(MS(k,j));
    end
end

MZ = sqrt(N)*MZ/trace(MZ);


%*****************************************************************************************

X=zeros(M);

for s = 1:M
    for sp = 1:M
        temp_y = randn;
        temp_z = randn;
        X(s,sp) = (temp_y + 1.*imag1*temp_z);
    end
end


rho_temp=X*X';

norm1=zeros(1,M);

for j = 1:N
    for k = 1:N
        s=N*(j-1)+k;
        
        ttt=rho_temp(s,s)/MZ(j,k);
        
        for sp = 1:M
            X(sp,s)=X(sp,s)/sqrt(ttt);
        end
        
    end
end


rho_temp=X*X';
rho_temp = sqrt(M)*rho_temp/trace(rho_temp);


rho=rho_temp;
%rho = Decoh(rho_temp,M,p); % decoherence of the random state

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
    sx=N*(Evals_re(nn)+1)/1.414;
    sy=N*Evals_im(nn)/1.414;
    %sx=Evals_re(nn)+1;
    %sy=Evals_im(nn);
    fprintf(hull,'%12.8f %16.12f\n',sx,sy);
end


fclose(hull);


%quit() %raskomentirovat dlja klastera
