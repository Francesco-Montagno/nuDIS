import lhapdf

pdf_set_name = f"new_PDF4LHC21_300"

pset = lhapdf.getPDFSet(pdf_set_name)
pdfs = pset.mkPDFs()


p_central = pdfs[0]