# primer3plus change log
## 1.0.6
**2019-11-11T08:52:02.431498**
removed unnecessay dependency

 - removed keats from development dependencies


## 1.0.5
**2019-10-07T15:34:50.778683**
bug fixes

 - design settings now restored with each run


## 1.0.4
**2019-10-07T15:18:25.381341**
improved error messages




## 1.0.3
**2019-10-07T14:50:00.887961**
bug fixes around long primer use

 - fixes 'SEQUENCE' and 'OVERHANG' keys for long primers
 - fixes primer locations for long primers
 - fixes primer pair product size for long primers
 - fixes bug in which product size did not match expected size for long primers. Will now automatically readjust product size parameter for long primers


## 1.0.2
**2019-10-07T11:15:48.129415**
bug fix

 - primer3plus.utils.anneal now by default ignores the case of the strings


## 1.0.1
**2019-10-07T08:02:08.512973**
fix overhang bug

 - Now, if overhangs are already explicitly set, the overhangs will be concatented if long_primer and overhangs are being used


## 1.0.0
**2019-10-06T10:48:12.224597**
stable release




## 1.0.0b1
**2019-10-03T13:36:40.865734**
bug fix

 - fixed documentation and arguments for `presets.included`


## 1.0.0b0
**2019-10-03T12:06:44.931634**





## 1.0.0a1
**2019-10-03T11:48:16.918526**
py35, py36, py37 support




## 1.0.0a
**2019-10-03T10:21:58.167755**
initial release


