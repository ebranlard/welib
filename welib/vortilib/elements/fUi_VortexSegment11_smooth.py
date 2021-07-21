import numpy as np
    
def fUi_VortexSegment11_smooth(xa = None,ya = None,za = None,xb = None,yb = None,zb = None,visc_model = None,t = None,bComputeGrad = None): 
    # !!!!! No intensity!!!
    norm_a = np.sqrt(xa * xa + ya * ya + za * za)
    norm_b = np.sqrt(xb * xb + yb * yb + zb * zb)
    denominator = norm_a * norm_b * (norm_a * norm_b + xa * xb + ya * yb + za * zb)
    crossprod = np.array([[ya * zb - za * yb],[za * xb - xa * zb],[xa * yb - ya * xb]])
    if bComputeGrad:
        Uout = np.zeros((12,1))
    else:
        Uout = np.array([[0],[0],[0]])
    
    # check for singularity */
#      fprintf('denom norma normb, #f #f #f\n',denominator,norm_a,norm_b)
    if denominator < 1e-17 and (visc_model < 4 or visc_model > 5):
        print('exit1')
        return Uout
    else:
        if (norm_a < 1e-08 or norm_b < 1e-08) and (visc_model < 4):
            print('exit2')
            return Uout
        else:
            # viscous model */
            Kv = 1.0
            if 0 == (visc_model):
                Kv = 1.0
            else:
                if 1 == (visc_model):
                    norm_r0 = np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb))
                    h = np.sqrt(crossprod(1) * crossprod(1) + crossprod(2) * crossprod(2) + crossprod(3) * crossprod(3)) / norm_r0
                    if (h < t):
                        Kv = h * h / t / t
                    else:
                        Kv = 1.0
                else:
                    if 2 == (visc_model):
                        norm_r0 = np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb))
                        h = np.sqrt(crossprod(1) * crossprod(1) + crossprod(2) * crossprod(2) + crossprod(3) * crossprod(3)) / norm_r0
                        Kv = 1.0 - np.exp(- 1.25643 * h * h / t / t)
                    else:
                        if 3 == (visc_model):
                            norm_r0 = np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb))
                            h = np.sqrt(crossprod(1) * crossprod(1) + crossprod(2) * crossprod(2) + crossprod(3) * crossprod(3)) / norm_r0
                            #             h = (norm_a+norm_b)/2;
                            Kv = h * h / np.sqrt(t ** 4 + h ** 4)
                        else:
                            if 4 == (visc_model):
                                norm_r0 = np.sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb))
                                Kv = 1.0
                                # delta*norm(r0)^2 */
                                denominator = denominator + t * norm_r0
                            else:
                                if 5 == (visc_model):
                                    Kv = 1.0
                                    # (delta l_0)^2 */
                                    denominator = denominator + t
                                else:
                                    if 33 == (visc_model):
                                        #             norm_r0 = sqrt((xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb));
#             h = sqrt(crossprod(1)*crossprod(1)+crossprod(2)*crossprod(2)+crossprod(3)*crossprod(3))/ norm_r0; # orthogonal distance r1xr2/r0 */
                                        h = (norm_a + norm_b) / 2
                                        Kv = h * h / np.sqrt(t ** 4 + h ** 4)
            #      fprintf('Kv #f\n',Kv)
            Kv = Kv / 4.0 / pi * (norm_a + norm_b) / denominator
            #      fprintf('Kv #f',Kv)
#      fprintf('crossprod #f\n',crossprod)
            Uout[np.arange[1,3+1]] = Kv * crossprod
            #      fprintf('Uout #f\n',Kv*crossprod)
            if bComputeGrad:
                d = (norm_a * norm_b + xa * xb + ya * yb + za * zb)
                if np.abs(d) > 1e-09:
                    ra = np.array([[xa],[ya],[za]])
                    rb = np.array([[xb],[yb],[zb]])
                    D = - (ra / norm_a ** 3 + rb / norm_b ** 3) - (1 / norm_a + 1 / norm_b) * 1 / d * (ra / norm_a * norm_b + rb / norm_b * norm_a + ra + rb)
                    GradRaRb = np.zeros((3,3))
                    GradRaRb[1,2] = za - zb
                    GradRaRb[2,1] = zb - za
                    GradRaRb[1,3] = yb - ya
                    GradRaRb[3,1] = ya - yb
                    GradRaRb[2,3] = xa - xb
                    GradRaRb[3,2] = xb - xa
                    GradU = 1 / (4 * pi) * 1 / d * (- (1 / norm_a + 1 / norm_b) * GradRaRb + crossprod * (np.transpose(D)))
                    Uout[3 + [np.arange[1,3+1]]] = GradU(1,np.arange(1,3+1))
                    Uout[6 + [np.arange[1,3+1]]] = GradU(2,np.arange(1,3+1))
                    Uout[9 + [np.arange[1,3+1]]] = GradU(3,np.arange(1,3+1))
            if sum(np.isnan(Uout)) > 0:
                print('Error fUi_Vortexline')
                #kbd
    
    #    printf("#4.3f #4.3f #4.3f #4.3f #4.3f\n",Uout(1),Uout(2),Uout(3),Kv,denominator); */
    return Uout