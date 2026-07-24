# Focal-Mechanism-From-SDR2MT
Continued the calculation and `decomposition of moment tensors` based on the previous repository <br>
-----
Vavryčuk (2015) discusses the decomposition of a general seismic moment tensor into three elementary components: the isotropic component (ISO), the double-couple component (DC), and the compensated linear vector dipole component (CLVD). This decomposition is mainly used to classify seismic source types and to provide a physical interpretation of non-double-couple mechanisms.

- The isotropic component represents a purely volumetric source. It corresponds to an equal expansion or contraction in all directions.

- The double-couple component is the most familiar component in earthquake seismology. It represents shear slip on a planar fault in an isotropic medium.

- The CLVD component is a non-double-couple deviatoric component. Unlike ISO and DC, it does not correspond to a simple elementary physical source by itself. However, it is mathematically necessary to complete the decomposition of a general moment tensor.

In anisotropic media, however, the moment tensor and source tensor may have different eigenvectors and different decompositions. As a result, moment tensor decomposition alone may be misleading. A non-zero ISO or CLVD component in the moment tensor may reflect anisotropic elastic properties rather than true tensile opening, closing, or source complexity. Therefore, Vavryčuk argues that in anisotropic focal regions, physical interpretation should be based on source tensor decomposition rather than moment tensor decomposition.

### **Refer to**: <br>
>> Vavryčuk, V., 2015. Moment tensor decompositions revisited, Journal of Seismology, 19(1), 231-252, doi: 10.1007/s10950-014-9463-y. <br>
<img width="1350" height="900" alt="Fig1" src="https://github.com/user-attachments/assets/26e0f821-7a7e-4ded-bf4b-476e1a74d0d4" />
<img width="1050" height="1020" alt="Fig2" src="https://github.com/user-attachments/assets/364f6980-dd9d-47dc-9b8b-69899351b8ce" />
Points in the first and third quadrants, where ISO and CLVD have the same sign: Positive ISO and CLVD indicate tensile opening, while negative ISO and CLVD indicate compressive closing. <br>
<img width="1050" height="1020" alt="Fig5" src="https://github.com/user-attachments/assets/546c93b1-6dee-4171-a3c2-5b948dde42e6" />

The third part of the script is used to determine which type of seismic event a given moment tensor most likely belongs to, after you already have the moment tensor results. For example:  
---
- Explosion events (EXPL);  
- Collapse events (CLAP);  
- Earthquake (EQs)
### **Refer to**: <br>
>> Sean R Ford, Gordon D Kraft, Gene A Ichinose, Seismic moment tensor event screening, Geophysical Journal International, Volume 221, Issue 1, April 2020, Pages 77–88, https://doi.org/10.1093/gji/ggz578 <br>
>> Mark J Hoggard, Janice L Scealy, Brent G Delbridge, Seismic moment tensor classification using elliptical distribution functions on the hypersphere, Geophysical Journal International, Volume 237, Issue 1, April 2024, Pages 1–13, https://doi.org/10.1093/gji/ggae011 <br>
>> Alvizuri, C., Silwal, V., Krischer, L., & Tape, C. (2018). Estimation of full moment tensors, including uncertainties, for nuclear explosions, volcanic events, and earthquakes. Journal of Geophysical Research: Solid Earth, 123, 5099–5119. https://doi.org/10.1029/2017JB015325 <br>
<img width="1050" height="1020" alt="Fig4" src="https://github.com/user-attachments/assets/38d78bdc-51e7-4126-9bb0-9674865aea3d" />


