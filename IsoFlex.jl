module IsoFlex
using FFTW, Statistics
export flexural, viscoelastic_lithos, viscoelastic_mantle
let
    global flexural
function flexural(ero::Array{Float64,2}; Te=30e3, dy=100,dx=100, 
    E = 100e9, g = 9.81, v = 0.25, pm = 3300, 
    pc = 2750, Nx = 0, Ny = 0, Nxy = 0, buffer::Int64 = 0,ncores::Int64 = 1)
    pw = 00;# water density
    m,n = size(ero);
    D = E*(Te)^3/(12*(1-v^2));
    Ly = m * dy;
    Lx = n * dx;
    bufferx::Int64 = 0
    buffery::Int64 = 0
    FFTW.set_num_threads(ncores)
    if mod(n,2) == 1
        bufferx .= bufferx.+1;
    end
    if mod(m,2) == 1
        buffery .= buffery.+1;
    end
    
    appndl = zeros(m,buffer.+bufferx);
    for i = 1:buffer +bufferx
       appndl[:,i] .=  ero[:,1];
    end
    appndr = zeros(m,buffer);
    for i = 1:buffer 
        appndr[:,i] .=  ero[:,end];
    end
    ero = [appndl ero appndr];
    appndu = zeros(buffer+buffery,n+2*buffer +bufferx);
    for i = 1:buffer+buffery
        appndu[i,:] .=  ero[1,:];
    end
    appndd = zeros(buffer,n+2*buffer+bufferx);
    for i = 1:buffer
        appndd[i,:] .=  ero[end,:];
    end
    ero = [appndu; ero; appndd];
    #imshow(ero)
    m,n = size(ero);
    k = zeros(Int(m/2+1),Int(n/2+1));
    k[1,1] = pc/(pm-pw);#(1+(D/(g*(pm-pc))*(2*pi*(1)/(Ly)).^4));
    dGR = g*(pm-pw)
    print(Nxy)
    for i = 2:Int(ceil((m-1)/2)+1)
        ky = (i - 1) / Ly
        k[i, 1] =  pc / (pm - pw) * 1. / (1. + ( D / dGR * (2 * pi * ky) .^ 4) + 4. / (pm * g) * pi ^ 2. * Ny * ky ^ 2);
    end
    for j = 2:Int(ceil((n-1)/2) + 1)
        kx = (j-1) / Lx
        k[1, j] = pc / (pm - pw) * 1. / ( 1. +  ( D / dGR * (2 * pi * kx) .^ 4) + 4. / (pm * g) * pi ^ 2. * Nx * kx ^ 2);
    end
    for i = 2:Int(ceil((m-1)/2)+1)
        for j = 2:Int(ceil((n-1)/2)+1)
           
            ky = (i-1) / Ly
            kx = (j-1) / Lx
            k[i,j] = pc/(pm - pw) .* 1. ./ (1. +  D / dGR * ( 2. * pi * sqrt( ky^2 + kx ^ 2. )) .^ 4. + 4. * pi ^ 2 / 
                    (pm*g) * ( Nx * kx ^ 2 + Ny * ky ^ 2 + Nxy * ky * kx ) );
        end
    end

        k = hcat(k, reverse(k[:,2:end-1],dims = 2));
        k = vcat(k, reverse(k[2:end-1,:],dims = 1));
    h = fft(ero.+1e-6,[1,2]);
    
    w2 = k .* real(h)+imag(h)*1im.*k;
    w = real(ifft(w2,[1,2]));
    w3 = w[buffer+1:end-buffer-bufferx,buffer+1:end-buffer-buffery];
    return w3
end
#global emaps,isomaps
end


let emaps = zeros(1,1,1),isomaps=zeros(1,1,1)#This emulates static behavior in julia
    global viscoelastic_lithos
    function viscoelastic_lithos(ero; Te=30e3, t=0, tt = 20e6, dt = 1e6, t_c = 100e6, dy=1000,dx=100, 
        E = 100e9, g = 9.81, v = .25, pm = 3300, pc = 2750, buffer::Int64 = 0, ncores::Int64=1)
        w = flexural(ero;Te=Te,E=E,v=v,g=g,pm=pm,pc=pc,buffer=buffer,ncores=ncores,dy=dy,dx=dx)
        wx = flexural(ero;Te=1e3,E=E,v=v,g=g,pm=pm,pc=pc,buffer=buffer,ncores=ncores,dy=dy,dx=dx)
        println(mean(ero))
        

        m,n = size(ero)
        w2=zeros(size(ero))
        if t == 0
            emaps = zeros(m, n, Int(floor(tt/dt)+1))
        else
            emaps[:,:,Int32(t/dt)+1] = w
            
        end
        if t == 0
            isomaps = zeros(m, n, Int(floor(tt/dt)+1))
        else
            isomaps[:,:,Int32(t/dt+1)] = wx
        end

        for t_i = dt:dt:t

            ds = -(emaps[:,:, Int(t_i/dt)+1 ] - isomaps[:,:, Int(t_i/dt)+1 ]) .* (dt/t_c)

            w .+= ds
            w2 .+= ds
            emaps[:,:, Int(t_i/dt) .+ 1] .+= ds

        end
            
            #w = w[buffer+1:end-buffer-bufferx,buffer+1:end-buffer-buffery];
        return w, w2
        end
    end
end