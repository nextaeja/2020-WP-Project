% KE from expectation value of p^2. i.e. do *-k^2 rather than *1i*k in expectation p.
function KE = PrintKEfromp2()
    En = EnergyExpectation();
    fprintf('KE from <p^2> = %.16e', En);
    KE = En;
end