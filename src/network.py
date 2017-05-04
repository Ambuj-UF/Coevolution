# -*- coding: utf-8 -*-
# Copyright 2014 by Ambuj Kumar, Kimball-Braun lab group, University of Florida.
# All rights reserved. This code is part of the EvoIntNet distribution and governed
# by its license. Please see the LICENSE file that should have been included
# as part of this package.
#
# Bug reports welcome: ambuj@ufl.edu
#
# Creates Cystoscape loadable protein interaction network XML file


import time
import math
from matplotlib import pyplot
from numpy import arange





def _unlist(dataList):
    newList = list()
    for obj in dataList:
        for inObj in obj:
            newList.append(inObj)

    return newList


def _fill_count_id(totalNode):
    newTotalNode = dict()
    for i, node in enumerate(totalNode):
        newTotalNode[node] = (i)

    return newTotalNode


def _count_interations(correl, totalNode):
    int_pairs = _unlist([key.split("-") for key, val in correl.items()])
    interaction_count = dict()
    for node in totalNode:
        interaction_count[node] = (int_pairs.count(node))

    return interaction_count


def _edge_count_id(correl, totalNode):
    edge_id = dict()
    for i, (key, val) in enumerate(correl.items()):
        edge_id[key] = (len(totalNode) + 100 + i)

    return edge_id


def _degree_dist(interaction_count):
    ki_list = [val for key, val in interaction_count.items()]
    total_ki = sum(ki_list)
    new_ki_list = sorted([float(x)/total_ki for x in ki_list])
    return new_ki_list


def scatterplot_node(dist_data):
    y = [math.log(val) for i, val in enumerate(dist_data)][::-1]
    x = [i for i, val in enumerate(dist_data)]
    pyplot.plot(x,y,'b.')
    pyplot.xlim(min(x)-1,max(x)+1)
    pyplot.ylim(min(y)-1,max(y)+1)
    pyplot.xlabel('node')
    pyplot.ylabel('ln(k/ktotal)')
    pyplot.show()


def scatterplot_connection(frac_data):
    x = [math.log(key) for key, val in frac_data.items()]
    y = [math.log(val) for key, val in frac_data.items()]
    pyplot.plot(x,y,'b.')
    pyplot.xlim(min(x)-1,max(x)+1)
    pyplot.ylim(min(y)-1,max(y)+1)
    pyplot.xlabel('k')
    pyplot.ylabel('nk')
    pyplot.show()


def node_fraction(interaction_count):
    k_groups = set([val for key, val in interaction_count.items()])
    retdata = dict()
    for ks in k_groups:
        retdata[ks] = (float(len([key for key, val in interaction_count.items() if val == ks]))/len(interaction_count))

    return retdata



def write_network_xml(correl):
    with open("network.xml", "w") as fp:
        fp.write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n")
        fp.write("<graph id=\"308\" label=\"EvoIntNet PPI Network\" directed=\"0\" cy:documentVersion=\"3.0\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:cy=\"http://www.cytoscape.org\" xmlns=\"http://www.cs.rpi.edu/XGMML\">\n")
        fp.write("\t<att name=\"networkMetadata\">\n")
        fp.write("\t<rdf:RDF>\n")
        fp.write("\t<rdf:Description rdf:about=\"http://www.cytoscape.org/\">\n")
        fp.write("\t<dc:type>Protein-Protein Interaction</dc:type>\n")
        fp.write("\t<dc:description>N/A</dc:description>\n")
        fp.write("\t<dc:identifier>N/A</dc:identifier>\n")
        fp.write("\t<dc:date>%s %s</dc:date>\n" %(time.strftime("%d/%m/%Y"), time.strftime("%H:%M:%S")))
        fp.write("\t<dc:title>EvoIntNet PPI Network</dc:title>\n")
        fp.write("\t<dc:source>http://www.cytoscape.org/</dc:source>\n")
        fp.write("\t<dc:format>Cytoscape-XGMML</dc:format>\n")
        fp.write("\t</rdf:Description>\n")
        fp.write("\t</rdf:RDF>\n")
        fp.write("\t</att>\n")
        fp.write("\t<att name=\"shared name\" value=\"EvoIntNet PPI Network\" type=\"string\"/>\n")
        fp.write("\t<att name=\"name\" value=\"EvoIntNet PPI Network\" type=\"string\"/>\n")
        fp.write("\t<att name=\"selected\" value=\"1\" type=\"boolean\"/>\n")
        fp.write("\t<att name=\"__Annotations\" type=\"list\">\n")
        fp.write("\t<att name=\"__Annotations\" value=\"\" type=\"string\"/>\n")
        fp.write("\t</att>\n")
        fp.write("\t<att name=\"layoutAlgorithm\" value=\"Grid Layout\" type=\"string\" cy:hidden=\"1\"/>\n")
    
        totalNode = set(_unlist([key.split("-") for key, val in correl.items()]))
        nodeData = _fill_count_id(totalNode)
        interaction_count = _count_interations(correl, totalNode)
        
        for key, val in nodeData.items():
            fp.write("\t<node id=\"%s\" label=\"%s\">\n" %(val, key))
            fp.write("\t\t<att name=\"shared name\" value=\"%s\" type=\"string\"/>\n" %key)
            fp.write("\t\t<att name=\"name\" value=\"%s\" type=\"string\"/>\n" %key)
            fp.write("\t\t<att name=\"selected\" value=\"0\" type=\"boolean\"/>\n")
            fp.write("\t\t<att name=\"canonicalName\" value=\"%s\" type=\"string\"/>\n" %key)
            fp.write("\t\t<att name=\"species\" value=\"\" type=\"string\"/>\n")
            fp.write("\t\t<att name=\"ALS-Aliases\" type=\"list\">\n")
            fp.write("\t\t\t<att name=\"ALS-Aliases\" value=\"\" type=\"string\"/>\n")
            fp.write("\t\t</att>\n")
            fp.write("\t\t<att name=\"searchTerm\" value=\"notSearchTerm\" type=\"string\"/>\n")
            fp.write("\t\t<att name=\"nbrConnections\" value=\"%s\" type=\"integer\"/>\n" %interaction_count[key])
            fp.write("\t</node>\n")


        edge_id = _edge_count_id(correl, totalNode)

        for key, val in correl.items():
            fp.write("\t<edge id=\"%s\" label=\"%s\" source=\"%s\" target=\"%s\" cy:directed=\"0\">\n" %(edge_id[key], edge_id[key], nodeData[key.split("-")[0]], nodeData[key.split("-")[1]]))
            fp.write("\t\t<att name=\"shared name\" value=\"\" type=\"string\"/>\n")
            fp.write("\t\t<att name=\"shared interaction\" value=\"pp\" type=\"string\"/>\n")
            fp.write("\t\t<att name=\"name\" value=\"\" type=\"string\"/>\n")
            fp.write("\t\t<att name=\"selected\" value=\"0\" type=\"boolean\"/>\n")
            fp.write("\t\t<att name=\"interaction\" value=\"pp\" type=\"string\"/>\n")
            fp.write("\t\t<att name=\"TSI-Sentences\" type=\"list\">\n")
            fp.write("\t\t\t<att name=\"TSI-Sentences\" value=\"\" type=\"string\"/>\n")
            fp.write("\t\t</att>\n")
            fp.write("\t\t<att name=\"TSI-PubEntry\" type=\"list\">\n")
            fp.write("\t\t\t<att name=\"TSI-PubEntry\" value=\"\" type=\"string\"/>\n")
            fp.write("\t\t</att>\n")
            fp.write("\t</edge>\n")

        fp.write("</graph>\n")

    dist_data = _degree_dist(interaction_count)
    scatterplot_node(dist_data)
    frac_data = node_fraction(interaction_count)
    scatterplot_connection(frac_data)
    



def write_network_sif(correl):
    with open("network.sif", "w") as fp:
        for key, val in correl.items():
            fp.write("%s\tpp\t%s\n" %(key.split("-")[0], key.split("-")[1]))




