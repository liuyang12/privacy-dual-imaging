function [H]=generate_gradient_operator(N1,N2,gradient_order)
% generate gradient operator matrix
if nargin>2
    go=gradient_order;
else
    go=2; % default
end

h1=[-1,1];
h2=[-1;1];

% H1
H1=zeros(N1*(N2-1),N1*N2);
for xx=1:N1
    for yy=1:N2-1
        H1(xx+(yy-1)*N1,xx+(yy-1)*N1)=h1(1);
        H1(xx+(yy-1)*N1,xx+yy*N1)=h1(2);
    end
end
H1=sparse(H1);

% H2
H2=zeros((N1-1)*N2,N1*N2);
for xx=1:N1-1
    for yy=1:N2
        H2(xx+(yy-1)*(N1-1),xx+(yy-1)*N1)=h2(1);
        H2(xx+(yy-1)*(N1-1),xx+1+(yy-1)*N1)=h2(2);
    end
end
H2=sparse(H2);
H=[H1;H2];

if go==2
    h3=[1,-2,1];
    h4=[1;-2;1];
    h5=[-1,1;1,-1];
    % H3
    H3=zeros(N1*(N2-2),N1*N2);
    for x=1:N1
        for y=1:N2-2
            H3(x+(y-1)*N1,x+(y-1)*N1)=h3(1);
            H3(x+(y-1)*N1,x+y*N1)=h3(2);
            H3(x+(y-1)*N1,x+(y+1)*N1)=h3(3);
        end
    end
    H3=sparse(H3);

    % H4
    H4=zeros((N1-2)*N2,N1*N2);
    for x=1:N1-2
        for y=1:N2
            H4(x+(y-1)*(N1-2),x+(y-1)*N1)=h4(1);
            H4(x+(y-1)*(N1-2),x+1+(y-1)*N1)=h4(2);
            H4(x+(y-1)*(N1-2),x+2+(y-1)*N1)=h4(3);
        end
    end
    H4=sparse(H4);

    % H5
    H5=zeros((N1-1)*(N2-1),N1*N2);
    for x=1:N1-1
        for y=1:N2-1
            H5(x+(y-1)*(N1-1),x+(y-1)*N1)=h5(1);
            H5(x+(y-1)*(N1-1),x+1+(y-1)*N1)=h5(2);
            H5(x+(y-1)*(N1-1),x+y*N1)=h5(3);
            H5(x+(y-1)*(N1-1),x+1+y*N1)=h5(4);
        end
    end
    H5=sparse(H5);
    H=[H;H3;H4;H5];
end

end