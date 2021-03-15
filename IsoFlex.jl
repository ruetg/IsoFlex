module IsoFlex
using PyPlot

using FFTW
export flexural, viscoelastic_lithos, viscoelastic_mantle
let
    global flexural
function flexural(ero::Array{Float64,2}; Te=30e3, dy=100,dx=100, 
    E = 100e9, g = 9.81, v = .25, pm = 3300, 
    pc = 2750,buffer::Int64 = 0,ncores::Int64 = 1)
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
    for i = 2:Int(ceil((m-1)/2)+1)
          k[i,1] =  ((pc)/(pm-pw))*1/(1+(D/(g*(pm-pw))*(2*pi*(i-1)/(Ly)).^4));
    end
    for j = 2:Int(ceil((n-1)/2) + 1)
          k[1,j] = (((pc)/(pm-pw))*1/(1+(D/(g*(pm-pw))*(2*pi*(j-1)/(Lx)).^4)));
    end
    for i = 2:Int(ceil((m-1)/2)+1)
        for j = 2:Int(ceil((n-1)/2)+1)
            k[i,j] = (((pc)/(pm-pw)).*1/(1+(D/(g*(pm-pw))*(2*pi*sqrt(((i-1)/(Ly))^2+((j-1)/(Lx))^2)).^4)));
        end
    end

        k = hcat(k, reverse(k[:,2:end-1],dims = 2));
        k = vcat(k, reverse(k[2:end-1,:],dims = 1));
    println("here")
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
    function viscoelastic_lithos(ero; Te=30e3, t=0, tt = 100e3, dt = 1e3, t_c = 10e6, dy=100,dx=100, 
        E = 100e9, g = 9.81, v = .25, pm = 3300, 
        bufferx = 0, buffery = 0, pc = 2750, buffer = 0, ncores=1)
        pw = 00;# water density
        m,n = size(ero);
        D = E*(Te)^3/(12*(1-v^2));
        Ly = m * dy;
        Lx = n * dx;
        println(t)


        FFTW.set_num_threads(ncores)
        if mod(n,2) == 1
            bufferx = bufferx+1;
        end
        if mod(m,2) == 1
            buffery = buffery+1;
        end

        appndl = zeros(m,buffer+bufferx);
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
        println(ero[1,:])
        h=h = fft(ero);

        m,n = size(ero);
        k = zeros(Int(m/2+1),Int(n/2+1));

        k[1,1] = pc/(pm-pw);#(1+(D/(g*(pm-pc))*(2*pi*(1)/(Ly)).^4));
        for i = 2:Int(ceil((m-1)/2)+1)
              k[i,1] =  ((pc)/(pm-pw))*1/(1+(D/(g*(pm-pw))*(2*pi*(i-1)/(Ly)).^4));
        end
        for j = 2:Int(ceil((n-1)/2) + 1)
              k[1,j] = (((pc)/(pm-pw))*1/(1+(D/(g*(pm-pw))*(2*pi*(j-1)/(Lx)).^4)));
        end
        for i = 2:Int(ceil((m-1)/2)+1)
            for j = 2:Int(ceil((n-1)/2)+1)
                k[i,j] = (((pc)/(pm-pw)).*1/(1+(D/(g*(pm-pw))*(2*pi*sqrt(((i-1)/(Ly))^2+((j-1)/(Lx))^2)).^4)));
            end
        end

        k = hcat(k, reverse(k[:,2:end-1],dims = 2));
        k = vcat(k, reverse(k[2:end-1,:],dims = 1));
        for i = 1:t
        end
        w = k .* h;
        w = real(ifft(w));
        if t == 0
            emaps = zeros(m, n, Int(floor(tt/dt)+1))
        else
            emaps[:,:,Int32(t/dt+1)] = w
        end
        if t == 0
            isomaps = zeros(m, n, Int(floor(tt/dt)+1))
        else
            isomaps[:,:,Int32(t/dt+1)] = ero .* pc/pm
        end
        w[:] .= 0.0
        emaps[:, :, Int(t/dt+1)] = w
        isomaps[:, :, Int(t/dt+1)] = w .* pc/pm

        for t_i = dt:dt:t
            ds = (emaps[:,:, Int(t_i/dt+1) ] - isomaps[:,:, Int(t_i/dt+1) ]) .* (exp(t_i/t_c ) - exp((t_i-dt)/t_c))
            w .+= ds
            emaps[:,:, Int(t_i/dt+1) ] .+= ds
            
        end
            
            w = w[buffer+1:end-buffer-bufferx,buffer+1:end-buffer-buffery];
        return w
        end
    end
end