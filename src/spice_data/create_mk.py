import os
def create_mk():
    file_path =os.path.dirname(__file__)
    tls_file = "latest_leapseconds.tls"
    bsp_file = "de432s.bsp"
    mk_file = "ss_kernel.mk"

    if not os.path.exists(os.path.join(file_path,mk_file)):
        with open(os.path.join(file_path,mk_file),'x') as f:
            kernel1 = os.path.join(file_path,tls_file)
            kernel2 = os.path.join(file_path,bsp_file)

            if len(kernel1) > 50:
                kernel1 = kernel1[0:50]+"+'\n'"+kernel1[50:]
            if len(kernel2) > 50:
                kernel2 = kernel2[0:50]+"+'\n'"+kernel2[50:]

            f.write("\\begindata\n\nKERNELS_TO_LOAD=(\n'"+kernel1+"',\n'"+kernel2+"'\n)\n\n\\begintext")
