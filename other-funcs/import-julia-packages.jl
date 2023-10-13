import Pkg
allpacks = split(ARGS[1],",")
for pack in allpacks
	try
		Pkg.add(pack)
		println("Added $pack")
	catch
		println("Couldn't add $pack")
	end
end
