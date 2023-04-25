import os
def create_mk():
    file_path =os.path.dirname(__file__)
    tls_file = "latest_leapseconds.tls"
    bsp_file = "de432s.bsp"
    mk_file = "ss_kernel.mk"

    if not os.path.exists(os.path.join(file_path,mk_file)):
        with open(os.path.join(file_path,mk_file),'x') as f:
            f.write("\\begindata\n\nKERNELS_TO_LOAD=(\n'"+os.path.join(file_path,tls_file)+"',\n'"+os.path.join(file_path,bsp_file)+"'\n)\n\n\\begintext")
