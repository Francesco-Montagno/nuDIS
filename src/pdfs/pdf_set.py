import lhapdf

pdf_set_name = f"PDF4LHC21_40"

pset = lhapdf.getPDFSet(pdf_set_name)
pdfs = pset.mkPDFs()


p_central = pdfs[0]