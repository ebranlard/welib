import numpy as np
class KalmanModel():
    def __init__(self,StateModel='nt1_nx7',Qgen_LSS=True):
        self.StateModel=StateModel
        self.Qgen_LSS=Qgen_LSS
        self.ThrustHack=False

        if Qgen_LSS:
            self.ColMap={
              ' ut1    ' : ' TTDspFA_[m]                   ' ,
              ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
              ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
              ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
              ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
              ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
              ' Qgen   ' : ' 97*{GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
              ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
              ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
              ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
            }
        else:
            self.ColMap={
              ' ut1    ' : ' TTDspFA_[m]                   ' ,
              ' psi    ' : ' {Azimuth_[deg]} * np.pi/180   ' , # [deg] -> [rad]
              ' ut1dot ' : ' NcIMUTVxs_[m/s]               ' ,
              ' omega  ' : ' {RotSpeed_[rpm]} * 2*np.pi/60 ' , # [rpm] -> [rad/s]
              ' Thrust ' : ' RtAeroFxh_[N]                 ' ,
              ' Qaero  ' : ' RtAeroMxh_[N-m]               ' ,
              ' Qgen   ' : ' {GenTq_[kN-m]}  *1000         ' , # [kNm] -> [Nm]
              ' WS     ' : ' RtVAvgxh_[m/s]                ' ,
              ' pitch  ' : ' {BldPitch1_[deg]} * np.pi/180 ' , # [deg]->[rad]
              ' TTacc  ' : ' NcIMUTAxs_[m/s^2]             ' 
            }

        if self.StateModel=='nt1_nx8':
            self.sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sAug        = np.array(['Thrust','Qaero','Qgen','WS'])
            self.sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            self.sInp        = np.array(['pitch'])
            self.sStor       = np.array(['WS'])
            self.bWSInStates     = True
            self.bThrustInStates = True
        elif self.StateModel=='nt1_nx7':
            self.sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sAug        = np.array(['Thrust','Qaero','Qgen'])
            self.sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            self.sInp        = np.array(['pitch'])
            self.sStor       = np.array(['WS'])
            self.bWSInStates     = False
            self.bThrustInStates = True
        elif self.StateModel=='nt1_nx6':
            self.sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sAug        = np.array(['Thrust','Qaero'])
            self.sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            self.sInp        = np.array(['Qgen','pitch'])
            self.sStor       = np.array(['WS'])
            self.bWSInStates     = False
            self.bThrustInStates = True
        elif self.StateModel=='nt1_nx5':
            self.sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            self.sAug        = np.array(['Qaero'])
            self.sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            self.sInp        = np.array(['Thrust','Qgen','pitch'])
            self.sStor       = np.array(['Thrust','WS'])
            self.bWSInStates     = False
            self.bThrustInStates = False











