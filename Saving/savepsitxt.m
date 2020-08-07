function savepsitxt()
    global psi nx ny nz
    f=fopen("psi.txt", "wt");
    for x= 1:nx
       for y = 1: ny
          for z= 1:nz
             fprintf(f,"%e\n",real(psi(x,y,z)));
             fprintf(f,"%e\n",imag(psi(x,y,z)));
          end
       end
    end
    fclose(f);
end