from asyncio.windows_events import NULL
from pyvis.network import Network




# Build cell communication network with results


def generateNetwork(CombinedSamples,CommunicatedCollection,CommunicationMap=NULL,celltype=NULL,img=NULL):
    net = Network()

    for samplename,sampleinfo in CombinedSamples.items():
        net.add_node(samplename,label=samplename)
    for Communicatedid,Collections in CommunicatedCollection.items():
        receivers = Collections['Receivers']
        secreters = Collections['Secreters']

        for secreter in secreters:
            for receiver in receivers:
                net.add_edge(secreter, receiver,title = Communicatedid )

    net.show("test.html")



    





























