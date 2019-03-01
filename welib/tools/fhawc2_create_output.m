function fhawc2_create_output()
    Sigs= {
  'aero alfa 1',
  'aero inflow_angle 1',
  'aero vrel 1',
  'aero induc 4 1 1',
  'aero induc 4 1 2',
  'aero induc 3 1 1',
  'aero induc 3 1 2',
  'aero secforce 1 1',
  'aero secforce 1 2',
  'aero cl 1',
  'aero cd 1',
  'aero windspeed 4 1 1',
  'aero windspeed 4 1 2',
  'aero velocity 4 1 1',
  'aero velocity 4 1 2',
  'aero gamma 1'};

r_ref=[0.562 0.788 1.350 1.845 2.070];
vr=linspace(0.05,2.25,10);
vr=unique(sort([vr r_ref]));

vr=[
    0.1125
    0.3375
    0.5625
    0.7875
    1.0125
    1.2375
    1.4625
    1.6875
    1.9125
    2.1375];



for is=1:length(Sigs)
    for ir=1:length(vr)
        fprintf('%s %.3f ;\n',Sigs{is},vr(ir));
    end
end
