function F=Fext(t,args)
if(t>=0)
    F=args.F0*cos(args.omega*t);
else
    F=0;
end

end