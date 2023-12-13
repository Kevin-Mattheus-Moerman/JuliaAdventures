using XML

# Define file name
filenameOut = "sample2.xml"

# Define document level
doc = XML.Document()

# Define declaration: <?xml version="1.0" encoding="UTF-8"?>
declarationNode = XML.Declaration(version = "1.0", encoding = "UTF-8")
push!(doc, declarationNode)

# Define food node
foodNode = XML.Element("food")
push!(doc, foodNode)

# Define fruit node
fruitNode = XML.Element("fruit")
push!(foodNode, fruitNode)

# Add banana nodes
for q=1:5
    idStr = string(q)
    priceStr = string(rand())

    bananaNode = XML.Element("banana", XML.Text("Delicious"); id = idStr, state = "ripe", price = priceStr)
    bananaNode.attributes["id"]=idStr
    bananaNode.attributes["price"]=priceStr

    push!(fruitNode, bananaNode)
end

# Write to XML file
XML.write(filenameOut, doc)
