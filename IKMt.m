function [IKMT,z1,x1] = IKMt(frqs,J)
close all;
%
% Projected Kirchhoff migration for a homogeneous 2D waveguide with 
% Dirichlet/Dirichlet boundary conditions.
% 
% The method was introduced in [1] for the case where the array spans the
% whole depth of the waveguide.
% 
% [1] C. Tsogka, D. A. Mitsoudis, and S. Papadimitropoulos. 
% Selective imaging of extended reflectors in two-dimensional waveguides. 
% SIAM J. Imaging Sci., 6(4):2714{2739, 2013.
% 
% Call as :
%
% [IKMT,z1,x1] = IKMt(frqs,J)
%
%
% frqs: the frequency (or frequencies) we use in the computation
%       (enter an array for multiple frequencies).
%
% J   : selective imaging using subspace projection on the Jth singular 
%       vector.
%       Enter 0 or leave blank for no projection,
%       Enter a vector if you want to project on more than one singular 
%       vector.
% 
% Examples of use:
%       [IKMT,z1,x1] = IKMt(73); Uses a single frequency and no projection.
%       [IKMt,z1,x2] = IKMt(31.875:3.75:118.125); Uses multiple frequencies
%                      and no projection.
%       [IKMt,z1,x1] = IKMt(73,1); Uses a single frequency and projection
%                       on the first singular vector.
%       [IKMt,z1,x1] = IKMt(73,[1,3,5]); Uses a single frequency and 
%                      projection on the first, third and fifth singular 
%                      vectors (the result of each projection is summed.


tt= tic;

% Argument checks.
if nargin == 0
    fprintf(['No arguments inserted. Running with default options '...
             '(f=73 Hz and J=0).\n'])
    frqs = 73;
    J = 0;
elseif nargin == 1
	J = 0;
end

% Frequency check.
if(length(frqs) == 1 && mod(frqs,3.75) ==0)
    error(['Frequency selected is a mode cut-off frequency. '... 
            'Choose again.'])
elseif length(frqs) > 1
%   Remove any mode cut-off frequencies from the array    
    frqs = frqs(mod(frqs,3.75)~=0); 
end


% Waveguide characteristics.
D = 200;
c0 = 1500;

% Array specs.
h = 5;
xs(:,2) = 0 : h : D;
xs(:,1) = 0;
xr(:,2) = 0 : h : D;
xr(:,1) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Location of the scatterer. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0 = 380;
x0 = 100;
% Size of the scatterer (if not a point).
b = 40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatterer shape (options shown below). %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sc = 'V';

if upper(sc) =='P'
%   Point
    zsc = z0;
    xsc = x0;
elseif upper(sc) =='S'
%   Semicircle
    thstep = pi/32;
    theta = pi/2 : thstep : 3*pi/2;
    zsc = z0 + b/2*cos(theta);
    xsc = x0 + b/2*sin(theta);
elseif upper(sc) =='L'
%   L-shape
    ho = b/32;
    zsc = [z0-b/2*ones(size(0:ho:b-ho)), z0-b/2:ho:z0+b/2];
    xsc = [x0-b/2:ho:x0+b/2, x0+b/2*ones(size(0:ho:b-ho))];
elseif upper(sc) =='H'
%   Horizontal Screen
    ho = b/32;
    zsc = z0-b/2:ho:z0+b/2;
    xsc = x0*ones(size(0:ho:b));
elseif upper(sc) == 'V'
%   Vertical Screen
    ho = b/32;
    zsc = z0*ones(size(0:ho/2:b));
    xsc = x0-b/2:ho/2:x0+b/2;
end

if xsc(1) < 0 || xsc(end) > D
    error('Scatterer out of bounds. Redefine cross-range position.')
end

% Reference frequency and wavelength.
f0 = 75;
lambda0 = c0/f0;

% Imaging/Search Domain.
zmin = z0 - 5*lambda0;
zmax = z0 + 5*lambda0;
xmin = 0;
xmax = D;
z1 = zmin : 2 : zmax; 
x1 = xmin : 2 : xmax; 
[Z, X] = meshgrid(z1,x1);

% Initializing Matrices
IKMT = zeros(length(x1),length(z1));
Nw = length(frqs);
npm = floor((2 * D * frqs(end))/c0);
phat_proj = zeros(Nw,npm,npm);

% Matrix creating start
tm = tic;

% Create and project the response matrix for the object.
count = 1;
for frq = frqs
	kn = (2 * pi * frq)/c0;
	NPM = floor((2 * D * frq)/c0);
	        
%   Create the response matrix.    
    p = 0;
    for i = 1 : length(zsc)
        p = p + green_inf(xs(:,1),xs(:,2),zsc(i),xsc(i),kn,D,NPM)...
        *green_inf(zsc(i),xsc(i),xr(:,1),xr(:,2),kn,D,NPM).';
    end
    
%   Create matrix VV containing the vertical eigenvunfctions.
    m = 1:NPM;
    lambdam = (m*pi / D).^2;
    betam = sqrt(kn*kn - lambdam);
    VV = sqrt(2/D) *sin(xr(:,2)*sqrt(lambdam));
    
%   Creating the multiplying matrices.
    Dbinv = diag(betam);
    B = h*(VV'*VV);
    [Ua,Sa,~] = svd(B);
    Dvinv = diag(1./diag(Sa));
    SJ = Dvinv*Ua'*VV';
    Shat  = SJ*p*SJ';
    
%   Projection of the response matrix.
    ptemp = Dbinv*(Ua*Shat*Ua')*Dbinv;

    phat_proj(count,1:NPM,1:NPM) = ptemp;
    count = count + 1; 
end
% Matrix creating finish.
tM = toc(tm);

% Imaging start.
ti = tic;

count = 1;
fcount = 1;
for frq = frqs
%   Wavenumber.
	kn = (2 * pi * frq)/c0;
%   Number of Propagating Modes.
	NPM = floor((2 * D * frq)/c0);
	p = squeeze(phat_proj(fcount,1:NPM,1:NPM));
	fcount = fcount + 1;

%   Creating a filtered version of the projected matrix by 
%   using the subspace projection method.
    if J ~= 0
    %J : filter
        filteredpr = zeros(size(p));
        [U,S,V] = svd(p);
        for m = J
            filteredpr = filteredpr + U(:,m)*S(m,m)*V(:,m)';
        end
        p = filteredpr;
    end

%   Computing IKMT.
    for m = 1 : NPM
        lambdam = (m * pi) / D;
        betam = sqrt(kn*kn-lambdam*lambdam);
        Xm = @(xs) sqrt(2/D)*sin(m*pi*xs/ D);
        for n = 1 : NPM
           lambdan = (n * pi) / D;
           betan = sqrt(kn*kn-lambdan*lambdan);
           Xn = @(xr) sqrt(2/D)*sin(n*pi*xr/D);
           
           IKMT = IKMT + exp(-1i*(betam+betan)*abs(Z-xs(n,1)))...
                 .*Xn(X).*Xm(X)*p(m,n);
        end
    end
   
    if(length(frqs) > 1)
        fprintf('Run %d out of %d \n',count,Nw)
    end
    count = count + 1;
end

% Imaging finish.
tI = toc(ti);

% Plot the modulus of the image.
figure
imagesc(z1,x1,abs(IKMT));
axis equal
axis tight
colormap jet;
hold on
if(upper(sc) =='P')
    plot(zsc,xsc,'w*')
else
    plot(zsc,xsc,'w')
end
drawnow;
shg;

% Autosave a .fig and .eps file of the image.
if upper(sc) == 'P'

    if length(frqs)==1
        filename = ['IKMt_P_f',num2str(frq),'_zsc_'...
                    ,num2str(z0),'_xsc_',num2str(x0),'_inf.fig'];
        saveas(gcf,filename)
    
        filename = ['IKMt_P_f',num2str(frq),'_zsc_'...
                 ,num2str(z0),'_xsc_',num2str(x0),'_inf.eps'];
        saveas(gcf,filename,'psc2')
    else
        filename = ['IKMt_P_f',num2str(frqs(1)),...
                    '_to_',num2str(frqs(end)),'_zsc_',num2str(z0),'_xsc_'...
                    ,num2str(x0),'_inf.fig'];
        saveas(gcf,filename)
    
        filename = ['IKMt_P_f',num2str(frqs(1)),...
                    '_to_',num2str(frqs(end)),'_zsc_',num2str(z0),'_xsc_'...
                    ,num2str(x0),'_inf.eps'];
        saveas(gcf,filename,'psc2')
    end
    
else

    if length(frqs)==1
        filename = ['IKMt_',upper(sc),num2str(b),'_f',num2str(frq),...
                    '_zsc_',num2str(z0),'_xsc_',num2str(x0),'_inf.fig'];
        saveas(gcf,filename)
    
        filename = ['IKMt_',upper(sc),num2str(b),'_f',num2str(frq),...
                    '_zsc_',num2str(z0),'_xsc_',num2str(x0),'_inf.eps'];
        saveas(gcf,filename,'psc2')
    else
        filename = ['IKMt_',upper(sc),num2str(b),'_f',num2str(frqs(1)),...
                    '_to_',num2str(frqs(end)),'_zsc_',num2str(z0),...
                    '_xsc_',num2str(x0),'_inf.fig'];
        saveas(gcf,filename)
    
        filename = ['IKMt_',upper(sc),num2str(b),'_f',num2str(frqs(1)),...
                    '_to_',num2str(frqs(end)),'_zsc_',num2str(z0),...
                    '_xsc_',num2str(x0),'_inf.eps'];
        saveas(gcf,filename,'psc2')
    end
end


tT = toc(tt);

% Print computation times.
fprintf('\n\n       Excecution time breakdown    \n')
fprintf('--------------------------------------- \n')
fprintf('Total time                   : %7.5f \n',tT)
fprintf('Creating Data and Projecting : %7.5f \n',tM)
fprintf('Image computation            : %7.5f \n',tI)
fprintf('Miscellaneous                : %7.5f \n', tT-tM-tI)


end
