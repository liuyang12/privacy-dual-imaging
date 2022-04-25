function [x]=ccgi(H,P,y,x_0,mu_0,rho,tol,min_iter,plot_flag,N1,N2,pinv_mat)
% min ||Hx||_1
% s.t.Px=y;

%%
% x=x_0;
x=pinv(P)*y;
w=H*x;
r=zeros(size(w));
s=zeros(size(y));
object_kp=norm(w(:),1);
o=object_kp;
mu=mu_0;
mu_bar=1e10;
converged=false;
iter=0;
max_iter=500;
% pinv_mat=pinv(H'*H+P'*P);
if plot_flag==1
figure;
norm_w_hx=norm(w-H*x);
norm_px_y=norm(P*x-y);
norm_x=[norm(x)];
norm_w=[norm(w)];
norm_r=[norm(r)];
norm_s=[norm(s)];
end

while ~converged
    iter=iter+1
    
    % update w
    tmp=H*x-1/mu*r;
    flag1=tmp>1/mu;
    flag2=tmp<-1/mu;
    w=(tmp-1/mu).*flag1+(tmp+1/mu).*flag2;
    
    % update x
    x=pinv_mat*(H'*(w+1/mu*r)+P'*(y-1/mu*s));
    
    % update r & s
    r=r+mu*(w-H*x); 
    s=s+mu*(P*x-y);
     
    % update mu
    mu=min(mu_bar,mu*rho);   
    
    % stop criterion
    object_k=object_kp;
    object_kp=norm(H*x,1);
    stop_c1=object_k-object_kp;
    stop_c2=norm(object_k-object_kp)/norm(object_kp);
    % [stop_c1,stop_c2];
    if ((stop_c1<tol && stop_c1>=0) || abs(stop_c1)<tol/100)   && stop_c2<tol  && iter>min_iter
        converged=true;
        disp('iteration is converged');
    elseif iter>max_iter
        converged=true;
        disp('iter number reached maximun number');
    end    
    
    if plot_flag==1
%             norm_x=[norm_x,norm(x)];
%             norm_w=[norm_w,norm(w)];
%             norm_r=[norm_r,norm(r)];
%             norm_s=[norm_s,norm(s)];
%             norm_w_hx=[norm_w_hx norm(w-H*x)]
%             norm_px_y=[norm_px_y norm(P*x-y)]
            o=[o,norm(H*x,1)]
%             subplot(241);plot(norm_w_hx);title('w-hx');
%             subplot(242);plot(norm_px_y);title('px-y');
%             subplot(243);plot(norm_x);title('x');
%             subplot(244);plot(norm_w);title(['w' sprintf(' mu %1.3f',mu)]);
%             subplot(245);plot(norm_r);title('norm r');
%             subplot(246);plot(norm_s);title('norm s');
            subplot(121);plot(o);title('object value');
            subplot(122);imshow(reshape(x,[N1,N2]));title('rst img');
            pause(0.01);
    end
end
if plot_flag==1
    close;
end
end