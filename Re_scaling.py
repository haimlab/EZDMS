def ReScaledEnrichmentRatio(SiReported,SMedianSynonymous,SMedianBottom):
    SiScaled = ((SiReported - SMedianSynonymous)/-SMedianBottom)+1
    return SiScaled

def calculateEnrichmentRatio(frX,frWT,urX,urWT):
    EnrichmentRatio = ((frX/frWT)/(urX/urWT))
    return EnrichmentRatio

def calculatePreference(EnrichmentRatioList,PreferenceSite):
    Preference = EnrichmentRatioList[PreferenceSite]/sum(EnrichmentRatioList)
    return Preference
