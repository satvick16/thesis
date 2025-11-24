coxn = 0.01973
coxp = 0.018665

u0n = 0.04359
u0p = 0.00432

vtn = 0.471
vtp = 0.423

vgs = 0.75
vdd = 0.75

length = 45e-9

# driver
rdrv = 56

wln_driver = 1 / (rdrv * u0n * coxn * (vgs - vtn))
wlp_driver = 1 / (rdrv * u0p * coxp * (vgs - vtp))

wn_driver = wln_driver * length
wp_driver = wlp_driver * length

# replica
rrep = 280

wln_replica = 1 / (rrep * u0n * coxn * (vgs - vtn))
wlp_replica = 1 / (rrep * u0p * coxp * (vgs - vtp))

wn_replica = wln_replica * length
wp_replica = wlp_replica * length

print("Driver NMOS width (m):", wn_driver)
print("Driver PMOS width (m):", wp_driver)
print("Replica NMOS width (m):", wn_replica)
print("Replica PMOS width (m):", wp_replica)
