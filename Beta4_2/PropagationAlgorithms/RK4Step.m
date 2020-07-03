% Alderwick's RK4 method
function psiStepped = RK4Step()
    global psi V dt;
    
    start_potential = V;
	halfway_potential = V;
	final_potential = V;
    
    grad1 = dbdt2(psi,start_potential);
	
	psiT = psi + grad1 * (dt / 2);
	grad2 = dbdt2(psiT,halfway_potential);
	
	psiT = psi + grad2 * (dt / 2);
	grad3 = dbdt2(psiT,halfway_potential);
	
	psiT = psi + grad3 * (dt);
	grad4 = dbdt2(psiT,final_potential);
	
	psiStepped = psi + dt / 6 * (grad1 + 2 * grad2 + 2 * grad3 + grad4);
    
end
function updatedPsi = dbdt2(psi, V)
    global kSquared mass hBar;
    
    fftPsi = fftn(psi);
    fftPsi = -1i*hBar/(2*mass)*kSquared.*fftPsi;
    updatedPsi = ifftn(fftPsi);
    updatedPsi = updatedPsi - 1i*V.*psi/hBar;
end