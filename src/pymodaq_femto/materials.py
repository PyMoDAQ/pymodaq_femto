from pypret.material import SellmeierF1



FS_extended = SellmeierF1(coefficients=[0.0000000, 0.6961663,
                               0.0684043, 0.4079426,
                               0.1162414, 0.8974794,
                               9.8961610],
                 freq_range=[1e-7, 6.7e-6],
                 name="FS",
                 long_name="Fused silica (fused quartz) extended range")