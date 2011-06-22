#!/usr/bin/env python

from database import DatabaseConnection


def main():
    """Main function"""
    connection = DatabaseConnection()
    connection.cursor()
    
    results = connection.get_chip_interactions()
    print results

if __name__ == "__main__":
    main()