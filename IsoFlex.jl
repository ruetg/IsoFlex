module IsoFlex

using FFTW
export flexural

function flexural(ero; Te=30e3,dy=100,dx=100, 
    E = 100e9, g = 9.81, v = .25, pm = 3300, 
    bufferx = 0, buffery = 0, pc = 2750,buffer = 0)
    pw = 00;# water density
    m,n = size(ero);
    D = E*(Te)^3/(12*(1-v^2));
    Ly = m * dy;
    Lx = n * dx;
    if mod(n,2) == 1
        bufferx = bufferx+1;
    end
    if mod(m,2) == 1
        buffery = buffery+1;
    end
    
    appndl = zeros(m,buffer+bufferx);
    for i = 1:buffer +bufferx
       appndl[:,i] =  ero[:,1];
    end
    appndr = zeros(m,buffer);
    for i = 1:buffer 
        appndr[:,i] =  ero[:,end];
    end
    ero = [appndl ero appndr];
    appndu = zeros(buffer+buffery,n+2*buffer +bufferx);
    for i = 1:buffer+buffery
        appndu[i,:] =  ero[1,:];
    end
    appndd = zeros(buffer,n+2*buffer+bufferx);
    for i = 1:buffer
        appndd[i,:] =  ero[end,:];
    end
    ero = [appndu; ero; appndd];
    h = fft(ero);
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
    k = [k  reverse(k[:,2:end-1],dims = 2)];
    k = [k; reverse(k[2:end-1,:],dims = 1)];
    w = k .* h;
    w = real(ifft(w));
    w = w[buffer+1:end-buffer-bufferx,buffer+1:end-buffer-buffery];
    return w
end
end