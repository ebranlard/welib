This folder contains examples of the "yams symbolic" framework which uses Kanes method to form the equations of motion.

There are three layers for the library:
- L1: low level layer, offer the most level of flexibility, but requires more lines of code
- L2: intermediate level, uses notions of RigidBodies, and flexibible bodies. The user has to specify the connection between the bodies, the degrees of freedom and the loads acting on them.
- L3: predefined models defined in welib/yams/models/ such as wind turbine models (FTNSB) or one rigid body models



# Examples of model names for L3 wind turbine models:
- Rigid body: F2T0RNA, F0T1N0S0, F6T0RNA 
- Flexible tower: F0T1RNA, F2T1RNA, F000000T2RNA,  etc

```python
model = get_model('F2T1RNA', mergeFndTwr=False, linRot=True,
                  yaw='zero', tilt='fixed',tiltShaft=True,
                  rot_elastic_type='SmallRot',
                  rot_elastic_smallAngle=False,
                  orderMM=1,
                  orderH=1,
                 )

model = get_model('F000000T2RNA', mergeFndTwr=True, linRot=False,
                  yaw='zero', tilt='fixed',tiltShaft=True,
                  rot_elastic_type='Body', #rot_elastic_type='Body', 'Body' or 'SmallRot'
                  rot_elastic_smallAngle=False, # Keep False nu needs to be second order
                  #rot_elastic_smallAngle=False,
                  orderMM=1,
                  orderH=1,
                  twrDOFDir=['x','y','x','y'], # Order in which the flexible DOF of the tower are set
                 )


#model = get_model_one_body('B100101', linRot=False, orderMM=1, orderH=1, )

```
