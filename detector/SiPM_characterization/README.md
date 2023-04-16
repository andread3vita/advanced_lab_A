# SiPM
## Description
The silicon photomultiplier (SiPM) is a radia- tion detector with extremely high sensitivity, high efficiency, and very low time jitter. It is based on reversed biased p/n diodes and can directly detect light from near ultra violet to near infrared. SiPMs are employed in all those applications where low light/radiation levels must be measured and quantified with high precision.

A SiPM consists of a matrix of small-sized sensitive elements called micro-cells (or pixels) all connected in parallel. Each micro-cell is a Geiger-Mode avalanche photo-diode (GM-APD) working beyond the breakdown voltage (Vbd) and it integrates a resistor for passive quenching.

In quiescent mode, the diode is reversed biased to Vbias = Vbd + Vov (Vov is the overvoltage, i.e. the excess bias beyond Vbd) and no current flows through the circuit. Once triggered, the avalanche process is self-sustaining meaning that, without quenching, a steady current flows indefinitely in the device (this is the main purpose of the quenching resistor Rq).

There are two main working phases: discharge phase and recovery phase. The first one correspons to the rising edge of the signal (time constant Cd x Rs, with Rs as the internal diode resistor), while the slower trailing edge is the recovery phase with time constant Cd x Rq. The amplitude of the SiPM pulse increases with the overvoltage while both the rising time and the recovery time are mainly determined by Cd, Rs, and Rq.

During recovery, the GM-APD cannot detect other photons. For this reason, a GM-APD cannot count the number of incoming photons unless the photon rate is lower than the inverse of the recharge time. SiPMs overcome this limita- tion thanks to the parallel arrangement of several micro-cells. When N photons are detected (which means that N photons arrive on N differ- ent micro-cells producing N single-cell signals) the SiPM output pulse is N-times larger than the single-cell response. In that, the N independent current pulses just add up at the SiPM termi- nals. Note that both the amplitude and the area of each SiPM pulse, which is the total charge delivered by the detector, are proportional to the number of detected photons (neglecting primary and correlated noise for the time being). If the number of incoming photons is larger than the number of micro-cells in the SiPM, saturation occurs and neither the amplitude nor the area of the output pulse can give information on the number of incoming photons anymore.

##  SiPM characterization

## References
* https://www.first-sensor.com/cms/upload/appnotes/AN_SiPM_Introduction_E.pdf
* https://www.ll.mit.edu/sites/default/files/publication/doc/geiger-mode-avalanche-photodiodes-three-aull-ja-7893.pdf
