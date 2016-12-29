def getCategory(line):
	if line.startswith('GPT'):
		return (255,0,0)
	elif line.startswith('PPT'):
		return (0,255,0)
	elif line.startswith('TPT'):
		return (0,0,255)
	elif line.startswith('XPT'):
		return (255,255,0)
	elif line.startswith('algae'):
		return (0,255,255)

	return (0,0,0)


if fontsize:
	monofont = ImageFont.truetype('FreeMono.ttf', fontsize)
	font_width, font_height = ImageDraw.Draw(Image.new('RGB', (100, 100), (0,0,0))).textsize('X', font=monofont)
	offset = (longest_acc + 1) * font_width
	font_height = int(font_height * 1.25)
	font_width = int(font_width * 1.4)
	height = len(order) * font_height
	width = len(next(iter(data.values()))) * font_width + offset

	im = Image.new('RGB', (width + 4*borders, height), (255,255,255))
	draw = ImageDraw.Draw(im)

	for y, elem in enumerate(order):
		if borders:
			draw.rectangle(((borders, y*font_height), (2*borders, (y+1)*font_height)), category[elem])
			draw.rectangle(((width + 2*borders, y*font_height), (width + 3*borders, (y+1)*font_height)), category[elem])
			if clusters:
				draw.rectangle(((0, y*font_height), (borders, (y+1)*font_height)), clusters[elem])
				draw.rectangle(((width + 3*borders, y*font_height), (width + 4*borders, (y+1)*font_height)), clusters[elem])
		draw.text((1 + 2*borders, y*font_height), elem, fill = (0,0,0), font = monofont)
		try:
			for x, value in enumerate(data[elem]):
				draw.rectangle(((offset + x*font_width + 2*borders, y*font_height), (offset + (x+1)*font_width + 2*borders, (y+1)*font_height)), value)
				draw.text((offset + x*font_width + 2*borders, y*font_height), alignments[elem][x], fill = (0,0,0), font = monofont)
		except KeyError:
			pass
else:
	height = len(order)
	width = len(next(iter(data.values())))

	im = Image.new('RGB', (width + 4*borders, height), (255,255,255))

	for y, elem in enumerate(order):
		for x in range(borders):
			im.putpixel((x + borders, y), category[elem])
			im.putpixel((x + width + 2*borders, y), category[elem])

			if clusters:
				im.putpixel((x, y), clusters[elem])
				im.putpixel((x + width + 3*borders, y), clusters[elem])
		try:
			for x, value in enumerate(data[elem]):
				im.putpixel((x + 2*borders, y), value)
		except KeyError:
			pass

if membranecoords:
	draw = ImageDraw.Draw(im, mode = 'RGBA')
	if fontsize:
		for x0, x1 in membranecoords:
			draw.rectangle(((offset + (x0-1)*font_width + 2*borders, 0), (offset + x1*font_width + 2*borders, height)), (0,0,0,120))
	else:
		for x0, x1 in membranecoords:
			draw.rectangle(((x0 + 2*borders, 0), (x1 + 2*borders + 1, height)), (0,0,0,120))

im.save(saveFN)
