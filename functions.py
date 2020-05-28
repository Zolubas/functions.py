import numpy as np
import cmath 
#-------------------------------------------------------#
#                                                       #
#               ASSOCIAÇÕES SIMPLES                     #
#                                                       #
#-------------------------------------------------------#
def divisor_de_tensao(tensao,z_quero,z_outra):
    return tensao*z_quero/(z_quero+z_outra)

def serie(z1,z2): #2 impedancias monofasicas em serie
    return z1 + z2
def paralelo(z1,z2): #2 impedancias monofasicas em paralelo
    return (z1*z2/(z1+z2))

#-------------------------------------------------------#
#                                                       #
#                   NUMEROS COMPELXOS                   #
#                                                       #
#-------------------------------------------------------#

def phase(x): #Rerturns phaase in degrees
    return ((180/np.pi)*(np.arctan2(x.imag,x.real)))

def modulo(x): #Returns module of complex number
    return np.sqrt(x.real*x.real + x.imag*x.imag)

def polar(x): #print phase/_angle of a fasor
    print(str(modulo(x)) + "/_" + str(phase(x)))

def polar_legenda(Nome_da_Variavel,x,unidade):
    print(str(Nome_da_Variavel) + " = " + str(modulo(x)) + "/_ " + str(phase(x))+ " " + str(unidade))

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

#-------------------------------------------------------#
#                                                       #
#                Trifásico e circuito CA                #
#                                                       #
#-------------------------------------------------------#
alpha = pol2cart(1,120)
def FP(complexo): #retorna fator potencia
    return complexo.real/np.sqrt(complexo.real**2 + complexo.imag**2)
def sequencia_direta(V):
    alpha = pol2cart(1,120)
    print("\nTensões em sequencia Direta\n")
    polar_legenda("V",V,"V")
    polar_legenda("V.alpha²",V*alpha**2,"V")
    polar_legenda("V.aplha",V*alpha,"V")
    return V,V*(alpha**2),V*alpha

def sequencia_inversa(V):
    alpha = pol2cart(1,120)
    print("\nTensões em sequencia Direta\n")
    polar_legenda("V",V,"V")
    polar_legenda("V.alpha²",V*alpha,"V")
    polar_legenda("V.aplha",V*alpha**2,"V")
    return V,V*(alpha),V*alpha**2

def sequencia_zero(V):
    return np.array([V,V,V])

def vlinha2vfase(VAB,VBC,VCA,sequencia_d_i):
    #Gerador simetrico e linha equiibrada
    if(sequencia_d_i == "d"):
        fator=pol2cart(np.sqrt(3),30)
        van = VAB/fator
        vbn = VBC/fator
        vcn = VCA/fator
        return van,vbn,vcn
    elif(sequencia_d_i == "i"):
        fator=pol2cart(np.sqrt(3),-30)
        van = VAB/fator
        vbn = VBC/fator
        vcn = VCA/fator
        return van,vbn,vcn
    else:
        print("\nVocê errou a sequência\n")

def vfase2vlinha(van,vbn,vcn,sequencia_d_i):
    if(sequencia_d_i == "d"):
        fator=pol2cart(np.sqrt(3),30)
        VAB = van*fator
        VBC = vbn*fator
        VCA = vcn*fator
        return VAB,VBC,VCA
    elif(sequencia_d_i == "i"):
        fator=pol2cart(np.sqrt(3),-30)
        VAB = van*fator
        VBC = vbn*fator
        VCA = vcn*fator
        return VAB,VBC,VCA
    else:
        print("\nVocê errou a sequência\n")

#-------------------------------------------------------#
#                                                       #
#           Sistemas de potência desequilibrados        #
#                                                       #
#-------------------------------------------------------#

def estrela_isolada_desequilibrada(van,vbn,vcn,zl,za,zb,zc):
    #Condições de USO
    #   Fonte simétrica equilibrada, Linha de transmissão equilibrada
    ya = 1/(zl + za)
    yb = 1/(zl + zb)
    yc = 1/(zl + zc)
    # o Correto é Vn'n
    VNNlinha = (van*ya + vbn*yb + vcn*yc)/(ya + yb + yc)

    #RESOLVENDO AS CORRENTES
    Vanbncn = np.array([van,vbn,vcn])
    vnnlinhavetor = np.array([VNNlinha,VNNlinha,VNNlinha])
    Y = np.matrix([[ya, 0, 0],[0, yb, 0],[0, 0, yc]])
    Iabc = np.dot(Y,Vanbncn- vnnlinhavetor) 
    Iabc_arr = [Iabc.item(0),Iabc.item(1),Iabc.item(2)]
    #Resolvendo as Tensões
    Z = np.matrix([[za, 0, 0],[0, zb, 0],[0, 0, zc]])
    Vanlbnlcnl = np.dot(Z,Iabc_arr) + vnnlinhavetor
    #Respostas
    print("\nRespostas Estrela isolada\n")
    polar_legenda("Ya ",ya,"siemens")
    polar_legenda("Yb ",yb,"siemens")
    polar_legenda("Yc ",yc,"siemens")
    polar_legenda("VN'N ",VNNlinha,"v")
    polar_legenda("Ia ",Iabc.item(0),"A")
    polar_legenda("Ib ",Iabc.item(1),"A")
    polar_legenda("Ic ",Iabc.item(2),"A")
    polar_legenda("Ia+Ib+Ic",Iabc.item(0)+Iabc.item(1)+Iabc.item(2),"A")
    print("Tensões de fase em relação ao nó n  ")
    polar_legenda("Va'n ",Vanlbnlcnl.item(0),"v")
    polar_legenda("Vb'n ",Vanlbnlcnl.item(1),"v")
    polar_legenda("Vc'n ",Vanlbnlcnl.item(2),"v")
    print("Tensões de fase realmente na carga em relação ao nó n' (N LINHA):  ")
    polar_legenda("Va'n' ",Vanlbnlcnl.item(0) - VNNlinha,"v")
    polar_legenda("Vb'n' ",Vanlbnlcnl.item(1) - VNNlinha,"v")
    polar_legenda("Vc'n' ",Vanlbnlcnl.item(2) - VNNlinha,"v")
    
def delta2estrelaABC(ZAB,ZBC,ZCA):
    #Retorna um vetor de impedâncias 
    ZA = ZAB*ZCA/(ZAB+ZBC+ZCA)
    ZB = ZBC*ZAB/(ZAB+ZBC+ZCA)
    ZC = ZBC*ZCA/(ZAB+ZBC+ZCA)
    print("\nImpedâncias estrela equivalente\n")
    polar_legenda("ZA ",ZA,"siemens")
    polar_legenda("ZB ",ZB,"siemens")
    polar_legenda("ZC ",ZC,"siemens")
    return ZA,ZB,ZC

def estrela2deltABBCCA(ZA,ZB,ZC):
    ZAB = (ZA*ZB + ZB*ZC + ZC*ZA)/ZC
    ZBC = (ZA*ZB + ZB*ZC + ZC*ZA)/ZA
    ZCA = (ZA*ZB + ZB*ZC + ZC*ZA)/ZB
    return ZAB,ZBC,ZCA


def estrela_aterrada_por_impedancia(van,vbn,vcn,zp,za,zb,zc,zn):
    #Condições de USO
    #   Fonte simétrica equilibrada, Linha de transmissão equilibrada
    
    #Resolvendo as correntes
    Vanbncn = np.array([van,vbn,vcn])
    Z = np.matrix([[zp + za + zn, zn, zn],[zn, zp + zb + zn, zn],[zn, zn, zp + zc + zn]])
    Iabc = np.dot(np.linalg.inv(Z),Vanbncn) 
    # print(type(Iabc)) , numpy.matrix
    In = Iabc.item(0) + Iabc.item(1) + Iabc.item(2)
    #Resolvendo VNN'
    VNNlinha = zn*In
    print("\nRespostas Estrela aterrada por impedância\n")
    polar_legenda("Ia ",Iabc.item(0),"A")
    polar_legenda("Ib ",Iabc.item(1),"A")
    polar_legenda("Ic ",Iabc.item(2),"A")
    polar_legenda("VNNlinha ",VNNlinha,"v")
    polar_legenda("Ineurtro",In,"A")

def estrela_aterrada_por_impedancia_mutua_e_groundImpedance(van,vbn,vcn,zp,za,zb,zc,zn,zm,zg,zmg):
    #Condições de USO
    #   Fonte simétrica equilibrada, Linha de transmissão equilibrada
    #   Francamente aterrado -> zn = 0
    #   Sem mutuas -> zm = 0
    #   Impedância de ground nula e mutuas de ground nulas -> zg=zmg=0

    #Impedância de terra não nula (ground)
    ZG = zg - 2*zmg
    #Resolvendo as correntes
    Vanbncn = np.array([van,vbn,vcn])
    ZGMG = np.matrix([[ZG,ZG,ZG],[ZG,ZG,ZG],[ZG,ZG,ZG]])
    Zm = np.matrix([[0,zm,zm],[zm,0,zm],[zm,zm,0]])
    Z = np.matrix([[zp + za + zn, zn, zn],[zn, zp + zb + zn, zn],[zn, zn, zp + zc + zn]]) + Zm +ZGMG
    Iabc = np.dot(np.linalg.inv(Z),Vanbncn) 
    # print(type(Iabc)) , numpy.matrix
    In = Iabc.item(0) + Iabc.item(1) + Iabc.item(2)
    #Resolvendo VNN'
    VNNlinha = zn*In
    print("\nRespostas Estrela aterrada por impedância\n")
    polar_legenda("Ia ",Iabc.item(0),"A")
    polar_legenda("Ib ",Iabc.item(1),"A")
    polar_legenda("Ic ",Iabc.item(2),"A")
    polar_legenda("VNNlinha ",VNNlinha,"v")
    polar_legenda("Ineurtro",In,"A")

#-------------------------------------------------------#
#                                                       #
#                        Pu                             #
#                                                       #
#-------------------------------------------------------#
def x_pu_real(x_pu,Sbase_do_sistema,Strafo):
    return x_pu*Sbase_do_sistema/Strafo

def trafo_monofrasico_em_pu(vn1_linha,vn2_linha,Strafo,XpercentualpuTRANSFORMADA):
    vn1 = vn1_linha
    vn2 = vn2_linha 
    #x = reatância em pu
    #Valroes de base, escolho a tensão de cada lado pra ser a base de cada lado
    #Trafo de potecnia, mesma potencia de base
    vbp = vn1
    vbs = vn2
    sb = Strafo
    #Bases deriadas
    zbp = vbp**2/sb
    zbs = vbs**2/sb
    ibp = sb/vbp
    ibs = sb/vbs
    v1 = vn1/vbp
    v2 = vn2/vbs
    print("\nTransformador em Pu\n")
    print("Tensoes de base")
    polar_legenda("Vbase_primario",vbp,"V")
    polar_legenda("Vbase_secudário",vbs,"V")
    print("Potência de base")
    polar_legenda("Sbase",sb,"VA")
    print("Impedâncias de base")
    polar_legenda("Zbase_primario",zbp,"ohms")
    polar_legenda("Zbase_secundário",zbs,"ohms")
    print("Correntes de base")
    polar_legenda("Ibase_primario",ibp,"A")
    polar_legenda("Ibase_secundario",ibs,"A")
    print("Tensoes no trafo em pu")
    polar_legenda("v1",v1,"pu")
    polar_legenda("v2",v2,"pu")
    print("Circuito equivalente em pu")
    print("\n-------------"+str(XpercentualpuTRANSFORMADA)+" pu ------------ \nv1                          v2\n------------------------------")

def pu_trifasico_carga_unica(Vb_linha,Sb_maior_pot_fornecida_para_carga):
    Ib = Sb_maior_pot_fornecida_para_carga/(Vb_linha*np.sqrt(3)) # (Sb3phi/3)/(VL/RAIZ(3))
    Zb = (Vb_linha**2)/Sb_maior_pot_fornecida_para_carga # (Vl/raiz(3))/(Sb3phi/3)
    print("\npu TRIFASICO CARGA UNICA\n")
    polar_legenda("Vbase",Vb_linha,"v")
    polar_legenda("Sbase",Sb_maior_pot_fornecida_para_carga,"VA")
    polar_legenda("Ibase",Ib,"A")
    polar_legenda("Zbase",Zb,"ohms")

def pu_trifasico_multi_carga(Vb_linha,Sb_potencia_da_fonte):
    Ib = Sb_potencia_da_fonte/(Vb_linha*np.sqrt(3)) # (Sb3phi/3)/(VL/RAIZ(3))
    Zb = (Vb_linha**2)/Sb_potencia_da_fonte # (Vl/raiz(3))/(Sb3phi/3)
    print("\npu TRIFASICO Multi carga\n")
    polar_legenda("Vbase",Vb_linha,"v")
    polar_legenda("Sbase",Sb_potencia_da_fonte,"VA")
    polar_legenda("Ibase",Ib,"A")
    polar_legenda("Zbase",Zb,"ohms")

def trafo_trifasico_delta_delta_pu(vn1_linha,vn2_linha,Strafo,XpercentualpuTRANSFORMADA):
    vn1 = vn1_linha
    vn2 = vn2_linha 
    #x = reatância em pu
    #Valroes de base, escolho a tensão de cada lado pra ser a base de cada lado
    #Trafo de potecnia, mesma potencia de base
    vbp = vn1
    vbs = vn2
    sb = Strafo
    #Bases deriadas
    zbp = vbp**2/sb
    zbs = vbs**2/sb
    ibp = sb/vbp
    ibs = sb/vbs
    print("\nTransformador em Pu\n")
    print("Tensoes de base")
    polar_legenda("Vbase_primario",vbp,"V")
    polar_legenda("Vbase_secudário",vbs,"V")
    print("Potência de base")
    polar_legenda("Sbase",sb,"VA")
    print("Impedâncias de base")
    polar_legenda("Zbase_primario",zbp,"ohms")
    polar_legenda("Zbase_secundário",zbs,"ohms")
    print("Correntes de base")
    polar_legenda("Ibase_primario",ibp,"A")
    polar_legenda("Ibase_secundario",ibs,"A")
