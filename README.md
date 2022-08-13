# Beam-element-stiffness-using-Matlab
function Ke=cotes_quadrature(K,aa,bb,N,nodes)
syms z;
mm=size(K);
yy=zeros(mm);
ss=zeros(mm);
dd=zeros(mm);
vv=zeros(mm);
ppp=zeros(mm);
xx=zeros(mm);
qq=zeros(mm);
endpts=zeros(mm);
pppp=zeros(mm);
g=zeros(mm);

% yy=0;dd=0;vv=0;xx=0;ss=0;qq=0;endpts=0;o=0
% N=(nodes-1)*ceil(N/(nodes-1))
N1=N+1;
x=linspace(aa,bb,N1)';
h=(x(N1,1)-x(1,1))/N;    
g=K;
NNN=(nodes-1)*N+1;
for ii=1:N
 ppp=subs(g,z,x(ii,1));
pppp=subs(g,z,x(ii+1,1));
endpts=endpts+ppp+pppp;
end
% if nodes==6
     for ii=2:5:NNN
            p =aa+(ii-1)*h/(nodes-1);
            yy=yy+subs(g,z,p);
          
    end
           for ii=5:5:NNN
            p =aa+(ii-1)*h/(nodes-1);
            ss=ss+subs(g,z,p);
           end
           for ii=3:5:NNN
            p =aa+(ii-1)*h/(nodes-1);
            dd=dd+subs(g,z,p);
           end
           for ii=4:5:NNN
            p =aa+(ii-1)*h/(nodes-1);
            vv=vv+subs(g,z,p);
           end
           for ii=6:5:NNN
            p =aa+(ii-1)*h/(nodes-1);
            xx=xx+subs(g,z,p);
        end

  Ke=vpa((5*(h/(nodes-1))/288)*(19*endpts+75*(yy+ss)+50*(dd+vv)));
% else
%    
% for ii=2:6:N
%             p =a+ii*h;
%             yy=yy+subs(g,z,p);
% end
%   for ii=6:6:N
%             p =a+ii*h;
%             ss=ss+subs(g,z,p);
%   end
%   for ii=4:6:N
%             p =a+ii*h;
%             dd=dd+subs(g,z,p);
%   end
%   for ii=3:6:N
%             p =a+ii*h;
%             vv=vv+subs(g,z,p);
%   end
%  for ii=5:6:N
%             p =a+ii*h;
%             xx=xx+subs(g,z,p);
%  end
%  for ii=7:6:N
%             p =a+ii*h;
%             qq=qq+subs(g,z,p);
%  end
%  Ke=vpa((6*h/840)*(41*endpts+216*(yy+ss)+27*(vv+xx)+272*dd))
% end
