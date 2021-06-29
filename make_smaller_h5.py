import lumicks.pylake as lk

file = lk.File('20201120-183136 Marker 1 DNA 1.h5')

file_out = 'chromatin_dna.h5'
file.save_as(file_out, omit_data={'Force HF/*', 'Photon count/*',
                                  'Trap position/*', 'Temperature/*',
                                  'Scan/*', 'Kymograph/*',
                                  'Force clamp/*', 'Diagnostics/*',
                                  'Confocal diagnostics/*'})