% KE from expectation value of p. Take this, and mod square to get p^2
function KE = PrintKEfromp()
    global mass;
    
    Mom = MomentumExpectation();
    KE = sum(Mom.^2)/(2*mass);
    
    fprintf('KE from <p>^2 = %.16e', KE);
end