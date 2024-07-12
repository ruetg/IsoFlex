## Example using Julia - see tests.ipynb

## Example using MATLAB (no viscoelstic support):

L = zeros(501,501);

L(200:300,200:300) = 1000; %Square load (m)

Te = 10e3; %effective elastic thickness of the lithosphere, m

dx=1000; %grid spacing, m

dy=1000; %grid spacing, m

Pc = 2700; % Crust density kg/m3

Buffer = 100; %Edge buffer size in pixels (to prevent artifacts on the edges - assume continuous load)

W = IsoFlex(L,Te,dx,dy,Pc,Buffer);

mesh(W)

