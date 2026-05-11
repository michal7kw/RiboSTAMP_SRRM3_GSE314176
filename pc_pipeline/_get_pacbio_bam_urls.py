#!/usr/bin/env python3
"""Query NCBI's SRA Data Locator (SDL) for the ORIGINAL PacBio BAM URLs
for the three GSE314176 Revio runs. The .sra archive that prefetch + sam-dump
produces is SRA-Normalized — it strips PacBio aux tags (qs, pw, ip, np, rq)
that skera + IsoQuant need. The original submitted BAM is on AWS S3 at a
different URL listed in SDL's response with supertype='Original'.

Outputs TSV: srr  filename  size_bytes  url
"""
import json
import sys
import urllib.request
import xml.etree.ElementTree as ET

SRRS = ["SRR36480452", "SRR36480453", "SRR36480454"]

print("srr\tfilename\tsize_bytes\turl")
for srr in SRRS:
    # Use the run_new XML endpoint which exposes 'Original' supertype + S3 URL
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc={srr}"
    with urllib.request.urlopen(url) as r:
        xml = r.read().decode()
    root = ET.fromstring(xml)
    for srafile in root.iter("SRAFile"):
        if srafile.attrib.get("supertype") == "Original":
            filename = srafile.attrib.get("filename")
            size = srafile.attrib.get("size")
            alts = srafile.find("Alternatives")
            s3_url = alts.attrib.get("url") if alts is not None else ""
            # Convert s3:// to https
            if s3_url.startswith("s3://"):
                bucket, key = s3_url[5:].split("/", 1)
                https = f"https://{bucket}.s3.amazonaws.com/{key}"
            else:
                https = s3_url
            print(f"{srr}\t{filename}\t{size}\t{https}")
