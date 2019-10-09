# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 13:13:09 2019

This code will search pubmed for your keyword and
then generate a full list of citations

result keys: 'Count', 'RetMax', 'RetStart', 
'IdList', 'TranslationSet', 'TranslationStack', 'QueryTranslation'

papers keys: 'PubmedArticle', 'PubmedBookArticle'

papers['PubmedArticle'][0]['MedlineCitation']['Article'].keys: 
    'ELocationID', 'ArticleDate', 'Language', 
    'Journal', 'ArticleTitle', 'Pagination', 'Abstract', 
    'AuthorList', 'GrantList', 'PublicationTypeList'
    
@author: nlama
"""

from Bio import Entrez
import json 

class PubMedQuery(object):
    
    def __init__(self,search):
        self.api_key = "1aacc2b4a9f4d9c730813db9a4da68a45009"
        self.email = 'nlama@email.unc.edu'
        results = self.search(search)
#        print(results.keys())
#        print(results['Count'])
        id_list = results['IdList']
        self.papers = self.fetch_details(id_list)
        
        ##Create Citation
  ###########Class Functions##########################
    
    def search(self, query):
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        handle = Entrez.esearch(db='pubmed', 
                                sort='relevance', 
                                retmax='20',
                                retmode='xml', 
                                term=query)
        results = Entrez.read(handle)
        return results
    
    def fetch_details(self, id_list):
        ids = ','.join(id_list)
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        results = Entrez.read(handle)
        return results
    
    def citationGenerator(self, papers):
        '''APA: Last, F. M. (Year Published). Article title. Journal Name,
           Volume(Issue), pp. Page(s). doi:# or Retrieved from URL
        '''
        self.citations = []
        
        for paper in papers['PubmedArticle']:
            try:
                authors = ""
                period = ". "
                L = ""
#                F = ""
                I = ""
                numAuthors = len(paper['MedlineCitation']['Article']['AuthorList'])
                for count, author in enumerate (paper['MedlineCitation']['Article']['AuthorList']):
                   L = author['LastName']
#                   F = author['ForeName']
                   I = author['Initials']
                   authors += "{0}, {1}.".format(L,I[0])
                   if count < numAuthors - 1:
                       authors += ", " #add comma if there are more authors
                try:
                    year = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
                except:
                    year = paper['MedlineCitation']['Article']['ArticleDate'][0]['Year']
                title = paper['MedlineCitation']['Article']['ArticleTitle']
                journal = paper['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
                try:
                    volume = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['Volume']
                except:
                    volume = ""
                    period = ""
                try:
                    issue = "({0})".format(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['Issue'])
                except:
                    issue = ""
                for item in paper['MedlineCitation']['Article']['ELocationID']:
                    if item.attributes['EIdType'] == 'doi':
                        doi = str(item)
                citation = "{0} ({1}). {2} {3}, {4}{5}{6}{7}".format(authors, year, title, journal, volume, issue, period, doi)
                self.citations.append(citation)
            except Exception as e:
                print("missing info:",e)
                print("for paper: ", paper['MedlineCitation']['Article']['ArticleTitle'], "\n")
                continue
           
            

if __name__ == '__main__':
    qObj = PubMedQuery("Kevin Weeks")
    papers = qObj.papers
    qObj.citationGenerator(papers)
    
for c in qObj.citations:
    print(c)
    print("\n")
#    for i, paper in enumerate(papers):
    # Pretty print the first paper in full to observe its structure
#print(json.dumps(papers['PubmedArticle'], indent=2, separators=(',', ':')))
