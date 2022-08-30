function csi_new=phaseOffset(csi,phase_diff_12,phase_diff_13)

amp = abs(csi);
phase = angle(csi);
phase_new = [phase(1,:); phase(2,:)-phase_diff_12; phase(3,:)-phase_diff_13];

csi_new = amp.*cos(phase_new)+amp.*1i.*sin(phase_new);

    
end