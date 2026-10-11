import zipfile,re,sys
import xml.etree.ElementTree as ET
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def read_xlsx(path):
    z=zipfile.ZipFile(path); ss=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',NS):
            ss.append(''.join(t.text or '' for t in si.iter('{%s}t'%NS['m'])))
    wb=ET.fromstring(z.read('xl/workbook.xml'))
    names=[s.get('name') for s in wb.find('m:sheets',NS)]
    out={}
    for i,n in enumerate(names,1):
        sh=ET.fromstring(z.read(f'xl/worksheets/sheet{i}.xml')); rows=[]
        for r in sh.iter('{%s}row'%NS['m']):
            row={}
            for c in r.findall('m:c',NS):
                ref=re.match(r'([A-Z]+)',c.get('r')).group(1); v=c.find('m:v',NS); t=c.get('t')
                if t=='inlineStr': val=''.join(x.text or '' for x in c.iter('{%s}t'%NS['m']))
                elif v is None: val=None
                elif t=='s': val=ss[int(v.text)]
                else: val=v.text
                row[ref]=val
            rows.append(row)
        out[n]=rows
    return out
