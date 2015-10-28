import re
from bs4 import BeautifulSoup
with open('./feature_table.html') as handle:
    soup = BeautifulSoup(handle, 'html.parser')

content = soup.find(attrs={'class': 'field-item'}).get_text().encode('utf-8')
appendicies = content[content.rindex('Feature keys reference'):]
feature_key_section = appendicies[0:appendicies.index('Summary of qualifiers')]
qualifiers_section = appendicies[appendicies.index('Summary of qualifiers'):]
qualifiers_section = qualifiers_section[0:qualifiers_section.index('Controlled vocabularies')]

def parse_data(fksl, section_key):
    fksd = {}
    current_fk = None
    current_key = None
    current_content = ""
    for line in fksl:
        if line.startswith(section_key):
            current_fk = line.strip().split()[-1]
            fksd[current_fk] = {}
            current_key = None
            current_content = ""

        if len(line.strip()) == 0:
            continue

        if current_fk is not None:
            if not line.startswith(' ') and not line.startswith('\t'):
                if line.startswith('Mandatory qualifiers'):
                    current_key, current_content = 'Mandatory qualifiers', line[len('Mandatory qualifiers'):]
                else:
                    try:
                        current_key, current_content = re.split('\s{3,}', line.strip())
                    except:
                        print line.strip()
                fksd[current_fk][current_key] = current_content
            else:
                fksd[current_fk][current_key] += ' ' + line.strip()
    return fksd


fksd = parse_data(feature_key_section.split('\n')[10:-1], 'Feature Key      ')
fqsd = parse_data(qualifiers_section.split('\n')[5:-1], 'Qualifier    ')

# post processing
for k1 in fksd:
    for k2 in fksd[k1]:
        if k2.lower() == 'optional qualifiers':
            kd2r = {}
            k2d = fksd[k1][k2]
            for line in k2d.split('/')[1:]:
                try:
                    (key, value) = line.split('=', 1)
                    value = value.strip()
                    if '"' == value[-1]:
                        value = value.strip('"')

                    kd2r[key] = value
                except:
                    kd2r[key] = None
            fksd[k1][k2] = kd2r


import json
with open('feature_table.json', 'w') as handle:
    json.dump(fksd, handle)
with open('qualifier_table.json', 'w') as handle:
    json.dump(fqsd, handle)
