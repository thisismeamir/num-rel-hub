def wrapper(text: str, wrapNumber: int):
    words = text.split()  # Split text into words
    wrapped = []  # List to hold the wrapped lines
    current_line = []  # Temporary list to hold words for the current line
    current_length = 0  # Track the length of the current line

    for word in words:
        # Check if adding this word exceeds the wrapNumber
        if current_length + len(word) + (1 if current_line else 0) <= wrapNumber:
            current_line.append(word)
            current_length += len(word) + (1 if current_line else 0)  # Account for the space
        else:
            # If it exceeds, finalize the current line and start a new one
            wrapped.append(" ".join(current_line))
            current_line = [word]  # Start the new line with the current word
            current_length = len(word)

    # Add the last line if it has words
    if current_line:
        wrapped.append(" ".join(current_line))

    return wrapped


