import numpy as np
from numpy.fft import fft2, ifft2
def kcalc(D,pc,pm,pw,g,m,n,Ly,Lx,Ny=0,Nx=0,Nxy=0):
    k = np.zeros([np.int32(m/2+1),np.int32(n/2+1)]);
    f = pc/(pm-pw) # inverse for ero
    k[0,0] = f;
    dGR = g*(pm-pw)
    for i in range(1, int(np.ceil((m-1)/2)+1)):
        ky = (i) / Ly
        k[i, 0] =  f * 1. / (1. + ( D / dGR * (2 * np.pi * ky) ** 4) + 4. / (pm * g) * np.pi ** 2. * Ny * ky ** 2);
    
    for j in range(1, int(np.ceil((n-1)/2) + 1)):
        kx = (j) / Lx
        k[0, j] = f * 1. / ( 1. +  ( D / dGR * (2 * np.pi * kx) ** 4) + 4. / (pm * g) * np.pi ** 2. * Nx * kx ** 2);
    
    for i in range(1,int(np.ceil((m-1)/2)+1)):
        for j in range(1,int(np.ceil((n-1)/2)+1)):
           
            ky = (i) / Ly
            kx = (j) / Lx
            k[i,j] = f * 1. / (1. +  D / dGR * ( 2. * np.pi * np.sqrt( ky ** 2 + kx ** 2. )) ** 4. + 4. * np.pi ** 2 / 
                    (pm*g) * ( Nx * kx ** 2 + Ny * ky ** 2 + Nxy * ky * kx ) );
    k = np.hstack([k, np.flip(k[:,1:-1],1)]);
    k = np.vstack([k, np.flip(k[1:-1,:],0)]);
    return k
def flexural(ero, Te=20e3, dy=1000,dx=1000, E = 100e9, g = 9.81, v = 0.25, pm = 3300, pc = 2750, Nx = 0, Ny = 0, Nxy = 0, buffer = 0):
    
    pw = 0;# water density
    m,n = np.shape(ero);
    s  = np.zeros(((m+2*buffer),(n+2*buffer)))
    if buffer>0:
        s[buffer:-buffer,buffer:-buffer] = ero
    else:
        s[:,:]=ero
    
    D = E*(Te)**3/(12*(1-v**2));
    
    Ly = (m + 2 * buffer) * dy - dy;
    Lx = (n + 2 * buffer) * dx - dx;
    
    
    m,n = np.shape(s)
    h = fft2(s + 1e-10)
    
    k = kcalc(D,pc,pm,pw,g,m,n,Ly,Lx,Ny=Ny,Nx=Nx,Nxy=Nxy)

    w2 = k * h
        
    w = np.real(ifft2(w2))
    if buffer>0:
        w=w[buffer:-buffer,buffer:-buffer]
    return w

def viscoelastic_lith(ero, Te=20e3, dy=1000,dx=1000,
                      E = 100e9, g = 9.81, v = 0.25, pm = 3300,
                      pc = 2750, Nx = 0, Ny = 0, Nxy = 0, buffer = 0,
                      T0 = 100e6, t=0, dt=1e6, def1=[], tdef=[] ):
    pw = 0  # water density
    m, n = np.shape(ero)
    s = np.zeros(((m + 2 * buffer), (n + 2 * buffer)))
    if buffer > 0:
        s[buffer:-buffer, buffer:-buffer] = ero
    else:
        s[:, :] = ero

    D = E * (Te) ** 3 / (12 * (1 - v ** 2))
    D2 = E * (2e3) ** 3 / (12 * (1 - v ** 2))

    Ly = (m + 2 * buffer) * dy - dy
    Lx = (n + 2 * buffer) * dx - dx

    m, n = np.shape(s)

    h = fft2(s + 1e-6)
    k = kcalc(D, pc, pm, pw, g, m, n, Ly, Lx, Ny=Ny, Nx=Nx, Nxy=Nxy)
    k2 = kcalc(D2, pc, pm, pw, g, m, n, Ly, Lx, Ny=Ny, Nx=Nx, Nxy=Nxy)
    w2 = k * h
    w2 = np.real(ifft2(w2))


    w_end = k2 * h
    w_end = np.real(ifft2(w_end))
    if buffer > 0:
        w_end = w_end[buffer:-buffer, buffer:-buffer]

    tdef.append(w_end.copy())
    def1.append(w2.copy())
    nts = int(T0 * 3 / dt)  # Minimal change after 3xT0
    if t / dt > nts:
        ts = int(t / dt - nts)
        tdef[np.floor(t / dt) - nts] = 0
        def1[np.floor(t / dt) - nts] = 0
    else:
        ts = 0
    for i in range(ts,int(t / dt)):
        w0 = -dt * (def1[i] - tdef[i]) / T0
        w2 += w0
        def1[i] += w0

    return w2, def1, tdef


def viscoelastic_mantle(ero, Te=20e3, dy=1000,dx=1000,
                      E = 100e9, g = 9.81, v = 0.25, pm = 3300,
                      pc = 2750, Nx = 0, Ny = 0, Nxy = 0, buffer = 0,
                      T0 = 5e3, t=0, dt=1e6, def1=[], tdef=[]):
    pw = 0  # water density
    m, n = np.shape(ero)
    s = np.zeros(((m + 2 * buffer), (n + 2 * buffer)))
    if buffer > 0:
        s[buffer:-buffer, buffer:-buffer] = ero
    else:
        s[:, :] = ero

    D = E * (Te) ** 3 / (12 * (1 - v ** 2))

    Ly = (m + 2 * buffer) * dy - dy;
    Lx = (n + 2 * buffer) * dx - dx;

    m, n = np.shape(s)

    h = fft2(s + 1e-6)
    k = kcalc(D, pc, pm, pw, g, m, n, Ly, Lx, Ny=Ny, Nx=Nx, Nxy=Nxy)
    w2 = k * h


    w2 = np.real(ifft2(w2))
    print(np.shape(w2))
    if buffer > 0:
        w2 = w2[buffer:-buffer, buffer:-buffer]
    print(t/dt)
    print(len(tdef))
    if len(tdef) <= t/dt:
        tdef.append(w2.copy())

        w2*=0 #We start at zero
        def1.append(0)
        nts = int(T0*5/dt) # Minimal change after 3xT0

        if t / dt > nts:
            ts = int(t / dt - nts)
            tdef[ts] = 0
            def1[ts] = 0
        else:
            ts = 0

        for i in range(ts,int(t / dt)):
            w0 = -dt * (def1[i] - tdef[i]) / T0
            w2 += w0
            def1[i] += w0
    else:
        tdef[-1] += w2 #If we're adding more deformation to the current timestep...
        w2 = None

    return w2, def1, tdef