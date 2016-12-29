use strict;
use Data::Dumper;
use Test::More;
use Config::Simple;
use Time::HiRes qw(time);
use Bio::KBase::AuthToken;
use Bio::KBase::workspace::Client;
use fba_tools::fba_toolsImpl;

# Some initialization upon loading
local $| = 1;
if (!defined($ENV{'KB_AUTH_TOKEN'})) {
	require "Bio/KBase/fbaModelServices/ScriptHelpers.pm";
	$ENV{'KB_AUTH_TOKEN'} = Bio::KBase::fbaModelServices::ScriptHelpers::getToken();
}
my $fba_client = FBAToolsClient->new($ENV{'KB_AUTH_TOKEN'},"https://kbase.us/services/ws");

# Outward facing function
sub call {
    my($function, $parameters) = @_;
    print "Function: ";
    print $function;
    print "\nParameters: ";
    my $output;
    my $key;
    foreach $key (keys $parameters) {
        print "$key is $parameters->{$key}\n";
    }
    $output =  $fba_client->run_function($function, $parameters);
    return %$output;
}

# The Client to the perl implementation fba_toolsImpl and fba_toolsServer
{
	package FBAToolsClient;
	use strict;
	use Test::More;
    sub new {
        my ($class,$token,$wsurl) = @_;
        my $object = fba_tools::fba_toolsImpl->new();
        my $self = {
            token => $token,
            ws_url => $wsurl,
            user_id => undef,
            ws_client => undef,
            obj => $object,
        };
        my $auth_token = Bio::KBase::AuthToken->new(token => $token, ignore_authrc => 1);
        $self->{user_id} = $auth_token->user_id();
        $self->{ws_client} = new Bio::KBase::workspace::Client($self->{ws_url},token => $self->{token});
        return bless $self, $class;
    }
    # make an API Call
    sub run_function {
		my($self,$function,$parameters) = @_; # $name,$tests,$fail_to_pass,$dependency) = @_;
        # Set up stuff for fba_toolsServer
		my $testctx = LocalCallContext->new($self->{token}, $self->{user},[{'service' => 'fba_tools', 'method' => $function, 'method_params' => [$parameters]}],$function);
		$fba_tools::fba_toolsServer::CallContext = $testctx;
        # Get Output
		my $output;
		# eval {
			if (defined($parameters)) {
				$output = $self->{obj}->$function($parameters);
			} else {
				$output = $self->{obj}->$function();
			}
		# };
		return $output;
	}
    sub workspace_url {
        # Use this function to change workspace systems
        my($self, $wsurl) = @_;
        $self->{ws_url} = $wsurl;
        $self->{ws_client} = new Bio::KBase::workspace::Client($self->{ws_url},token => $self->{token});
    }
}



# Package below is necessary for setting up fba_toolsServer properly
{
    package LocalCallContext;
    use strict;
    sub new {
        my($class,$token,$user,$provenance,$method) = @_;
        my $self = {
            token => $token,
            user_id => $user,
            provenance => $provenance,
            method => $method
        };
        return bless $self, $class;
    }
    sub user_id {
        my($self) = @_;
        return $self->{user_id};
    }
    sub token {
        my($self) = @_;
        return $self->{token};
    }
    sub provenance {
        my($self) = @_;
        return $self->{provenance};
    }
    sub method {
        my($self) = @_;
        return $self->{method};
    }
    sub authenticated {
        return 1;
    }
    sub log_debug {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
    sub log_info {
        my($self,$msg) = @_;
        print STDERR $msg."\n";
    }
}
1;