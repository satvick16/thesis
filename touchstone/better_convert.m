S = sparameters('C2M__Z100_IL14_WC_BOR_H_L_H_THRU.s4p');

[s_dd, s_dc, s_cd, s_cc] = s2smm(S.Parameters);

freq = S.Frequencies;

Sdiff = sparameters(s_dd, freq, 50);

rfwrite(Sdiff, 'C2M__Z100_IL14_WC_BOR_H_L_H_THRU_two_port.s2p');
