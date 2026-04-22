import lhapdf

pdfs      = None
p_central = None

def configure(pdf_set_name):
    global pdfs, p_central
    pset      = lhapdf.getPDFSet(pdf_set_name)
    pdfs      = pset.mkPDFs()
    p_central = pdfs[0]
