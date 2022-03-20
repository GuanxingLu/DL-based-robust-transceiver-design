function [H,u_ori,u,v,lambda,P,miu,A,beta,r] = robust_transceiver_beamforming(H,Gama,sigma2_noise,P0)
    

    K = size(H,3);
    N = size(H,2);
    M = size(H,1);
    
    H = permute(H,[2,1,3]);
    iter = 1;
    r = 0;

    u = randn(N,K) + 1i*randn(N,K);
    %compute u by heuristic
    for k = 1:K
        D = H(:,:,k)*H(:,:,k)';
        [eigvec,eigmat]=eig(D);
        eigenvalue = diag(eigmat);
        maxeig = max(eigenvalue);
        for i=1:N
            if maxeig==eigenvalue(i)
                break;
            end
        end
        u(:,k) =eigvec(:,i);
    end
    
    for k = 1:K
        u(:,k) = u(:,k)/norm(u(:,k));
    end
    u_ori = u;
    
    lambda = Gama*ones(K,1);
    for k = 1:K
        lambda(k) = lambda(k)/norm(H(:,:,k));
    end
    
    Fix_poit_iter_Num = 20;
    
    v = randn(M,K) + 1i*randn(M,K);
    f = 0;

    for k = 1:K
        f = f + norm(v(:,k))^2;
        v(:,k) = v(:,k)/norm(v(:,k));
    end

    mem = 0;

    t = 0;
    rem = [f];

    beta = r*ones(K,1);
    alpha = zeros(K,1);
    for i = 1:K
        alpha(i) = sqrt(norm(H(:,:,i)));
    end

    chg = 1;
    
    sigma2 = sigma2_noise*ones(K,1);

        mem = f;
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
        end
        
       r = (P0-ones(K,1)'/A*sigma2_noise*ones(K,1))/(ones(K,1)'/A*ones(K,1));
       r = abs(r);

        for i = 1:K
            B(i) = sigma2(i)*norm(u(:,i))*norm(u(:,i)) + r;
        end
       
        miu = A\B;

        for k = 1:K
            v(:,k) = sqrt(miu(k))*v(:,k);
        end

  
        f = sum(sqrt(real(miu)));
        rem = [rem,f];

    lambda = real(lambda);

    P = 0;

    for user_index = 1:K
                P = P + norm(v(:,user_index))^2;
    end
    
    disp(['whole_iter ' num2str(iter) ', r=' num2str(r)...
        ', dP=' num2str(abs(P-P0))]); 
            
        
        
    
        






