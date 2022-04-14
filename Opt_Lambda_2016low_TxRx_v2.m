function [H,u_ori,u,v,lambda,P,miu,A,r] = Opt_Lambda_2016low_TxRx_v2(H,Gama,sigma2,P0)
    

    K = size(H,3);%M*N*K
    N = size(H,2);
    M = size(H,1);
    
    H = permute(H,[2,1,3]);%N*M*K
    iter = 1;
    r = 0;

    lambda = ones(K,1);
    
    Fix_poit_iter_Num = 20;
    
    %v = ones(M,K) + 1i*ones(M,K);
    u = ones(N,K);
    u_ori=u;
    v = randn(M,K) + 1i*randn(M,K);

%     v(:,1) = [0.0701+0.7443i,-0.3386+0.0235i].';
%     v(:,2) = [-1.3709+2.0320i,-0.9692-0.1711i].';
    
%     u(:,1) = [-0.7423-0.1885i,-0.2951-0.5713i].';
%     u(:,2) = [0.7580-0.6429i,-0.1084+0.0209i].';

f = 0;

for k = 1:K
    f = f + norm(v(:,k))*norm(v(:,k));
    v(:,k) = v(:,k)/norm(v(:,k));
end

mem = 0;


t = 0;
rem = [f];

while(abs(f-mem)>0.0001)
%while(abs(f)>0.05)

%while(t < 1000)
    mem = f;
    t = t+1;
    for k = 1:K
        tmp = zeros(N,N);
        for j = 1:K
            if j ~= k
                tmp = tmp + H(:,:,k)*v(:,j)*v(:,j)'*H(:,:,k)';
            end  
        end 
        tmp = tmp + sigma2*eye(N);
        u_ori(:,k) = eye(N)/tmp*H(:,:,k)*v(:,k);
        %u(:,k) = u_ori(:,k)./abs(u_ori(:,k))./2;
        u(:,k) = u_ori(:,k)/norm(u_ori(:,k));
        %u(:,k) = u_ori(:,k)./abs(u_ori(:,k))./2
    end
    
  
    for n = 1:Fix_poit_iter_Num
        tmp = zeros(M,M);
        for j = 1:K
            tmp = tmp + lambda(j)*H(:,:,j)'*u(:,j)*u(:,j)'*H(:,:,j);
        end
        tmp = tmp + eye(M);
        Upsilon = eye(M)/tmp;
        for k = 1:K
            lambda(k) = (Gama/(1+Gama))/(u(:,k)'*H(:,:,k)*Upsilon*H(:,:,k)'*u(:,k));
        end
    end
    
    for k = 1:K
        tmp = zeros(M,M);
        for j = 1:K
            if j ~= k
                tmp = tmp + lambda(j)*H(:,:,j)'*u(:,j)*u(:,j)'*H(:,:,j);
            end
        end
        tmp = tmp + eye(M);
        tmp = eye(M)/tmp;
        
        v(:,k) = tmp*H(:,:,k)'*u(:,k);
        v(:,k) = v(:,k)/norm(v(:,k));
    end
    
    A = zeros(K,K);
    B = zeros(K,1);
      
    for i = 1:K
        for j = 1:K
            if j == i
                A(i,j) = (1/Gama)*(u(:,i)'*H(:,:,i)*v(:,i)*v(:,i)'*H(:,:,i)'*u(:,i));
            else
                A(i,j) = -(u(:,i)'*H(:,:,i)*v(:,j)*v(:,j)'*H(:,:,i)'*u(:,i));
            end
        end
        
        B(i,1) = sigma2*norm(u(:,i))*norm(u(:,i));
    end
    
    miu = A\B;
    
    for k = 1:K
        v(:,k) = sqrt(miu(k))*v(:,k);
    end
    
    f = sum(sqrt(real(miu)));
    rem = [rem,f];
end

lambda = real(lambda);

r = (P0-ones(K,1)'/A*sigma2*ones(K,1))/(ones(K,1)'/A*ones(K,1));
% if r<0
%     y=1;
% end
%r = abs(r);
r = real(r);

sigma2_noise = sigma2*ones(K,1);
for i = 1:K
    B(i) = sigma2_noise(i)*norm(u(:,i))^2 + r;
end

miu = A\B;

for k = 1:K
    v(:,k) = sqrt(miu(k))*v(:,k)/norm(v(:,k));
end

P = 0;
for user_index = 1:K
            P = P + norm(v(:,user_index))^2;
end

disp(['2016low_v2 r=' num2str(r) ', dP=' num2str(abs(P-P0))]); 
    
        






