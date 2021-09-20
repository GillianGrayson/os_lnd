%str1 = getenv('SGE_TASK_ID');
str1 = '1';
number = str2num(str1);


attrac2 = fopen('spectrum_10.dat','w');

N=40; % system size
M=N^2-1; % auxiliary size
MM=N*N;
imag1=sqrt(-1);

alpha=0.1;

Ns=400;
Hist=zeros(1,Ns);

%%% Step 1
%%% NB: here we complete the F-basis by a normalized identity matrix as
%%% F{1} - just for the future, but will not actually use it.

k=1;
F{1}=sparse(eye(N))/sqrt(N);
for i=1:N
    for j=i+1:N
        k=k+1;
        F{k}=sparse([i j],[j i],[1 1]/sqrt(2),N,N);
        k=k+1;
        F{k}=sparse([i j],[j i],-imag1*[1 -1]/sqrt(2),N,N);
    end
end

for i=1:N-1
    k=k+1;
    temp=zeros(1,i+1);
    temp(1:i)=ones(1,i);
    temp(i+1)=-i;
    F{k}=sparse([1:i+1],[1:i+1],temp/sqrt(i*(i+1)),N,N);
end

for rudya=1:100
    for oslob=1:10
        
        %step = 1
        
        %%% Step 2: Generate X-matrix
        X = zeros(M,M);
        
        %seed = number;
        %randn('state', seed);
        temp_y = randn(M,'double');
        
        %seed = 1001 + number;
        %randn('state', seed);
        temp_z = randn(M,'double');
        
        X = (temp_y + imag1*temp_z)/2.;
        
        %step = 2
        
        
        %%% Step 3^ calculate G-matrix
        
        G = zeros(M,M);
        G = X*X';
        G = N*G/trace(G);
        
        %Eval_G = eig(G); % eigenvalues of K-matrix
        
        %write eigenvalues of rho-matrix to txt-files
        
        %fname_G = strcat('eval_g_',str1,'.txt');
        %save(fname_G, 'Eval_G', '-ASCII', '-double')
        
        %step = 3
        
        
        
        
        %%% Step 3^ calculate GA-matrix
        
        k=1;
        GA=zeros(MM,MM);
        for ks=1:M
            for kt=1:M
                ZZ=zeros(MM,MM);
                ZZ=kron(transpose(F{ks}),F{kt});
                GA=GA+G(ks,kt)*ZZ;
            end
        end
        
        
        %%*******************************************NEW PART
        %******************************* generacia Hamiltonian
        
        for n=1:N
            for m=1:N
                y(n,m)=normrnd(0,1)+i*normrnd(0,1);
            end
        end
        
        yy=ctranspose(y);
        BS=(y+yy)/sqrt(2);
        %BS=sqrt(N/4)*BS/sqrt(trace(BS*BS'));
        BS=sqrt(100*N)*BS/sqrt(trace(BS*BS'));
        
        
        %***************************************
        AS=alpha*1i*BS;
        
        
        b=kron(II,AS)+kron(conj(AS),II);
        
        GA=GA-b;
        %%*************************************************** END OF NEW PART
        
        
        
        ZS = eig(GA);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        for ks=2:MM
            
            r=real(ZS(ks))*N;
            
            r1=imag(ZS(ks))*N;
            
            
            fprintf(attrac2,'%8.4f %16.12f\n',r,r1);
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %step = 4
        
        oslob
        
    end
    
    oslob
    rudya
    
    
    
    
    
end
%quit()
fclose(attrac2);
