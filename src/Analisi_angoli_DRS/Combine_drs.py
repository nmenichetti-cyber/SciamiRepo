import xml.etree.ElementTree as ET
import argparse

def combine_drs(drs1, drs2, output_file):
     tree1 = ET.parse(drs1) # Parse the first DRS XML file
     tree2 = ET.parse(drs2) # Parse the second DRS XML file
     root1 = tree1.getroot() # Get the root element of the first DRS XML.
     root2 = tree2.getroot() # Get the root element of the second DRS XML
     print(root1.tag) # Print the tag of the root element of the first DRS XML
     all_events_1 = root1.findall('Event') # Find all Event elements in the first DRS XML
     all_events_2 = root2.findall('Event') # Find all Event elements in the second DRS XML
     new_root = ET.Element('DRSOSC') # Create a new root element for the combined DRS XML
     new_root.extend(all_events_1) # Add all Event elements from the first DRS XML to the new root element
     for event in all_events_2: # Iterate through all Event elements in the second DRS XML
        event.find('Serial').text = str(int(event.find('Serial').text)+len(all_events_1))   # ora <Serial>vecchio+ultimo_1</Serial>
        new_root.append(event) # Append each Event element from the second DRS XML to the new root element
     new_tree = ET.ElementTree(new_root)
     ET.indent(new_tree, space='  ')     # Python >= 3.9, aggiunge indentazione
     new_tree.write(output_file, xml_declaration=True, encoding='ISO-8859-1') # Write the combined DRS XML to a new file with XML declaration and specified encoding


def split_drs(drs1,output_file1,output_file2):
       tree1 = ET.parse(drs1) # Parse the DRS XML file
       root1 = tree1.getroot() # Get the root element of the DRS XML
       all_events_1 = root1.findall('Event') # Find all Event elements in the DRS XML
       new_root = ET.Element('DRSOSC') # Create a new root element for the split DRS XML
       new_root2 = ET.Element('DRSOSC') # Create a second new root element for the second split DRS XML
       print(f"Numero totale eventi: {len(all_events_1)}")
       print(f"Half: {len(all_events_1) // 2}")
       print(f"Primo serial: {all_events_1[0].find('Serial').text}")
       print(f"Ultimo serial: {all_events_1[-1].find('Serial').text}")
       for event in all_events_1: # Iterate through all Event elements in the DRS XML
         if int(event.find('Serial').text) <= len(all_events_1) // 2:   # Se il Serial è minore o uguale a 1000, aggiungi l'evento al nuovo root
               new_root.append(event)
         else: # Altrimenti, aggiungi l'evento al secondo nuovo root
               new_root2.append(event)      
       new_tree = ET.ElementTree(new_root)
       new_tree2 = ET.ElementTree(new_root2)
       ET.indent(new_tree, space='  ')     # Python >= 3.9, aggiunge indentazione
       ET.indent(new_tree2, space='  ')     # Python >= 3.9, aggiunge indentazione
       new_tree.write(output_file1, xml_declaration=True, encoding='ISO-8859-1') # Write the split DRS XML to a new file with XML declaration and specified encoding
       new_tree2.write(output_file2, xml_declaration=True, encoding='ISO-8859-1') # Write the second split DRS XML to a new file with XML declaration and specified encoding


      
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge di file XML dalla DRS4 in un unico file ordinato cronologicamente.')
    parser.add_argument('files', nargs='+', help='Lista dei file XML da unire')
    parser.add_argument('-o', '--output',nargs=2, default='merged.xml', help='Nome del file di output (default: merged.xml)')
    args = parser.parse_args()

    split_drs(args.files[0], args.output[0], args.output[1])     