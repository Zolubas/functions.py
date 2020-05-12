def phase(x): #Rerturns phaase in degrees
    return ((180/np.pi)*(np.arctan2(x.imag,x.real)))

def modulo(x): #Returns module of complex number
    return np.sqrt(x.real*x.real + x.imag*x.imag)

def polar(x): #print phase/_angle of a fasor
    print(str(modulo(x)) + "/_" + str(phase(x)))

def polar_legenda(Nome_da_Variavel,x,unidade):
    print(str(Nome_da_Variavel) + " = " + str(modulo(x)) + "/_" + str(phase(x))+ " " + str(unidade))

def cart2pol(x, y): #transform cartesian notation to polar (fasor)
    rho = np.sqrt(x**2 + y**2)
    phi = 180*np.arctan2(y, x)/np.pi
    return(rho, phi)

def pol2cart(rho, phi): #transform polar notation in an complex number
    phi = em_rad(phi)
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x + y*1j

def em_rad(x): #Degrees to radianus
    return np.pi*x/180
