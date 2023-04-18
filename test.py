# Import modules
import tkinter as tk
import _sqlite3
import re

# Define Model class
class Model:
    def __init__(self, email):
        self.email = email
        # Create a database connection
        self.connection = _sqlite3.connect('email.db')

        # Create a cursor object
        self.cursor = self.connection.cursor()

        # Create a table called email with a column called address
        self.cursor.execute('CREATE TABLE IF NOT EXISTS email (email TEXT)')

    @property
    def email(self):
        return self.__email

    @email.setter
    def email(self, value):
        # Validate the email using a regular expression
        pattern = r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'
        if re.fullmatch(pattern, value):
            self.__email = value
        else:
            raise ValueError(f'Invalid email address: {value}')

    def save(self):
        # Save the email to the database
        self.cursor.execute('INSERT INTO email VALUES (?)', (self.email,))
        # Commit the changes to the database
        self.connection.commit()
    
    def close(self):
        # Close the connection and the cursor
        self.connection.close()
        self.cursor.close()

# Define View class
class View:
    def __init__(self, controller):
        # Create a root window
        self.root = tk.Tk()
        self.root.title('Email Validator')
        # Create an entry widget for entering an email
        self.entry = tk.Entry(self.root)
        self.entry.pack()
        # Create a label widget for showing messages
        self.label = tk.Label(self.root)
        self.label.pack()
        # Bind the return key to the controller's save method
        self.entry.bind('<Return>', controller.save)
        # Set a reference to the controller
        self.controller = controller

    def get_email(self):
        # Get the email from the entry widget
        return self.entry.get()

    def show_message(self, message, color):
        # Show a message in the label widget with a color
        self.label.config(text=message, fg=color)

    def clear_entry(self):
        # Clear the entry widget
        self.entry.delete(0, tk.END)

# Define Controller class
class Controller:
    def __init__(self, view, model):
        # Set references to the view and the model
        self.view = view
        self.model = model

    def save(self, event):
        # Get the email from the view
        email = self.view.get_email()
        try:
            # Set the email to the model
            self.model.email = email
            # Save the email using the model's method
            self.model.save()
            # Show a success message in the view
            self.view.show_message('Email saved successfully', 'green')
            # Clear the entry in the view
            self.view.clear_entry()
        except ValueError as e:
            # Show an error message in the view
            self.view.show_message(e, 'red')

# Create an instance of Model with an empty email
model = Model('')
# Create an instance of View with a reference to Controller
view = View(Controller)
# Create an instance of Controller with references to View and Model
controller = Controller(view, model)
# Run the app using View's root window
view.root.mainloop()