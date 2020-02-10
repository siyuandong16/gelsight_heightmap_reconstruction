function T = poly_nd(order,varargin)
% POLY_ND(order,x,y,z,...) returns the polynomial terms of N variables.
% The input can be column vectors of numbers, chars, or symbols.
% The terms are orgnized in ascending dimensions.i.e x then x,y then x,y,z.
% EXAMPLE 1
%   syms x,y,z,t
%   poly_nd(4,x,y)
% ans =
%   [ 1, x, x^2, x^3, x^4, y, x*y, x^2*y, x^3*y, y^2, x*y^2, x^2*y^2, y^3, x*y^3, y^4]

%   poly_nd(2,x,y,z)
% ans =
%     [ 1, x, x^2, y, x*y, y^2, z, x*z, y*z, z^2]
%
%   poly_nd(2,x,y,z,t)
% ans =
%     [ 1, x, x^2, y, x*y, y^2, z, x*z, y*z, z^2, t, t*x, t*y, t*z, t^2]

%     
% EXAMPLE 2
%   x=[0 1 2]'; %numbers vector
%   y=['y1'; 'y2'; 'y3']; % chars vector
%   z=sym('z',[1,3]); %symbols vector
%   poly_nd(2,x,y,z)
% ans = 
%     [ 1, 0, 0, y1,    0, y1^2, z1,    0, y1*z1, z1^2]
%     [ 1, 1, 1, y2,   y2, y2^2, z2,   z2, y2*z2, z2^2]
%     [ 1, 2, 4, y3, 2*y3, y3^2, z3, 2*z3, y3*z3, z3^2]
%
% see also: polyval_nd

varnum = length(varargin);

for k=1:varnum
    
    arg=varargin{k};
    X{k} = sym(['X__' num2str(k)]);
    
    if isnumeric(arg) || strcmp(class(arg),'sym')
        arg=arg(:);
        eval(['X__' num2str(k) '= arg;']);
        
    elseif ischar(arg)
        for j=1:size(arg,1)
            tmp(j,1) = sym(arg(j,:));
        end
        
        eval(['X__' num2str(k) '= tmp;']);
        

    else
        error('Input type of %s not support',class(arg));
    end
    row(k)=size(arg,1);
    
end
row = unique(row);
if numel(row) > 1
    error('length of the inputs should be the same');
end
I=ones(row,1);

% T=sort(poly_nd_helper(order,X{:}));
T=eval(poly_nd_helper(order,X{:}));


function M=poly_nd_helper(order,varargin)

if numel(varargin) == 0 || order == 0
    M=sym('I');
else
    M=[];
    X=varargin(1:end-1); %X1 -- X_n-1
    Xn = varargin{end}; %Xn
    for k=0:order
        M=[M, poly_nd_helper(order-k,X{:}).*Xn.^k ];%recursively build ND polynomial based on its N-1 dimension
    end
end
